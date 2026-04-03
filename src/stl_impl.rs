// Ported from https://www.netlib.org/a/stl
//
// Cleveland, R. B., Cleveland, W. S., McRae, J. E., & Terpenning, I. (1990).
// STL: A Seasonal-Trend Decomposition Procedure Based on Loess.
// Journal of Official Statistics, 6(1), 3-33.

#![allow(clippy::too_many_arguments)]

use super::float::Float;

fn pow2<T: Float>(x: T) -> T {
    x * x
}

fn pow3<T: Float>(x: T) -> T {
    x * x * x
}

pub fn stl<T: Float>(
    y: &[T],
    np: usize,
    ns: usize,
    nt: usize,
    nl: usize,
    isdeg: i32,
    itdeg: i32,
    ildeg: i32,
    nsjump: usize,
    ntjump: usize,
    nljump: usize,
    ni: usize,
    no: usize,
    rw: &mut [T],
    season: &mut [T],
    trend: &mut [T],
    work: &mut [T],
) {
    let n = y.len();
    let work_size = n + 2 * np;
    let (work1, work) = work.split_at_mut(work_size);
    let (work2, work) = work.split_at_mut(work_size);
    let (work3, work) = work.split_at_mut(work_size);
    let (work4, work5) = work.split_at_mut(work_size);

    let mut userw = false;
    let mut k = 0;

    loop {
        onestp(
            y, n, np, ns, nt, nl, isdeg, itdeg, ildeg, nsjump, ntjump, nljump, ni, userw, rw,
            season, trend, work1, work2, work3, work4, work5,
        );
        k += 1;
        if k > no {
            break;
        }
        for i in 0..n {
            work1[i] = trend[i] + season[i];
        }
        rwts(y, n, work1, rw);
        userw = true;
    }

    if no == 0 {
        for v in rw.iter_mut() {
            *v = T::one();
        }
    }
}

fn ess<T: Float>(
    y: &[T],
    n: usize,
    len: usize,
    ideg: i32,
    njump: usize,
    userw: bool,
    rw: &[T],
    ys: &mut [T],
    res: &mut [T],
) {
    if n < 2 {
        ys[0] = y[0];
        return;
    }

    let mut nleft = 0;
    let mut nright = 0;

    let newnj = njump.min(n - 1);
    if len >= n {
        nleft = 1;
        nright = n;
        let mut i = 1;
        while i <= n {
            let ok = est(
                y,
                n,
                len,
                ideg,
                T::from_usize(i),
                &mut ys[i - 1],
                nleft,
                nright,
                res,
                userw,
                rw,
            );
            if !ok {
                ys[i - 1] = y[i - 1];
            }
            i += newnj;
        }
    } else if newnj == 1 {
        // newnj equal to one, len less than n
        let nsh = len.div_ceil(2);
        nleft = 1;
        nright = len;
        for i in 1..=n {
            // fitted value at i
            if i > nsh && nright != n {
                nleft += 1;
                nright += 1;
            }
            let ok = est(
                y,
                n,
                len,
                ideg,
                T::from_usize(i),
                &mut ys[i - 1],
                nleft,
                nright,
                res,
                userw,
                rw,
            );
            if !ok {
                ys[i - 1] = y[i - 1];
            }
        }
    } else {
        // newnj greater than one, len less than n
        let nsh = len.div_ceil(2);
        let mut i = 1;
        while i <= n {
            // fitted value at i
            if i < nsh {
                nleft = 1;
                nright = len;
            } else if i > n - nsh {
                nleft = n - len + 1;
                nright = n;
            } else {
                nleft = i - nsh + 1;
                nright = len + i - nsh;
            }
            let ok = est(
                y,
                n,
                len,
                ideg,
                T::from_usize(i),
                &mut ys[i - 1],
                nleft,
                nright,
                res,
                userw,
                rw,
            );
            if !ok {
                ys[i - 1] = y[i - 1];
            }
            i += newnj;
        }
    }

    if newnj != 1 {
        let mut i = 1;
        while i <= n - newnj {
            let delta = (ys[i + newnj - 1] - ys[i - 1]) / T::from_usize(newnj);
            for j in i + 1..=i + newnj - 1 {
                ys[j - 1] = ys[i - 1] + delta * T::from_usize(j - i);
            }
            i += newnj;
        }
        let k = ((n - 1) / newnj) * newnj + 1;
        if k != n {
            let ok = est(
                y,
                n,
                len,
                ideg,
                T::from_usize(n),
                &mut ys[n - 1],
                nleft,
                nright,
                res,
                userw,
                rw,
            );
            if !ok {
                ys[n - 1] = y[n - 1];
            }
            if k != n - 1 {
                let delta = (ys[n - 1] - ys[k - 1]) / T::from_usize(n - k);
                for j in k + 1..=n - 1 {
                    ys[j - 1] = ys[k - 1] + delta * T::from_usize(j - k);
                }
            }
        }
    }
}

fn est<T: Float>(
    y: &[T],
    n: usize,
    len: usize,
    ideg: i32,
    xs: T,
    ys: &mut T,
    nleft: usize,
    nright: usize,
    w: &mut [T],
    userw: bool,
    rw: &[T],
) -> bool {
    let range = T::from_usize(n) - T::one();
    let mut h = (xs - T::from_usize(nleft)).max(T::from_usize(nright) - xs);

    if len > n {
        h += T::from_usize((len - n) / 2);
    }

    let h9 = T::from_f64(0.999) * h;
    let h1 = T::from_f64(0.001) * h;

    // compute weights
    let mut a = T::zero();
    for j in nleft..=nright {
        w[j - 1] = T::zero();
        let r = (T::from_usize(j) - xs).abs();
        if r <= h9 {
            if r <= h1 {
                w[j - 1] = T::one();
            } else {
                w[j - 1] = pow3(T::one() - pow3(r / h));
            }
            if userw {
                w[j - 1] *= rw[j - 1];
            }
            a += w[j - 1];
        }
    }

    if a <= T::zero() {
        false
    } else {
        // weighted least squares
        for j in nleft..=nright {
            // make sum of w(j) == 1
            w[j - 1] /= a;
        }

        if h > T::zero() && ideg > 0 {
            // use linear fit
            let mut a = T::zero();
            for j in nleft..=nright {
                // weighted center of x values
                a += w[j - 1] * T::from_usize(j);
            }
            let mut b = xs - a;
            let mut c = T::zero();
            for j in nleft..=nright {
                c += w[j - 1] * pow2(T::from_usize(j) - a);
            }
            if c.sqrt() > T::from_f64(0.001) * range {
                b /= c;

                // points are spread out enough to compute slope
                for j in nleft..=nright {
                    w[j - 1] *= b * (T::from_usize(j) - a) + T::one();
                }
            }
        }

        *ys = T::zero();
        for j in nleft..=nright {
            *ys += w[j - 1] * y[j - 1];
        }

        true
    }
}

fn fts<T: Float>(x: &[T], n: usize, np: usize, trend: &mut [T], work: &mut [T]) {
    ma(x, n, np, trend);
    ma(trend, n - np + 1, np, work);
    ma(work, n - 2 * np + 2, 3, trend);
}

fn ma<T: Float>(x: &[T], n: usize, len: usize, ave: &mut [T]) {
    let newn = n - len + 1;
    let flen = T::from_usize(len);

    // get the first average
    let mut v: T = x.iter().take(len).copied().sum();
    ave[0] = v / flen;

    if newn > 1 {
        for (m, aj) in ave.iter_mut().take(newn).skip(1).enumerate() {
            // window down the array
            v = v - x[m] + x[len + m];
            *aj = v / flen;
        }
    }
}

fn onestp<T: Float>(
    y: &[T],
    n: usize,
    np: usize,
    ns: usize,
    nt: usize,
    nl: usize,
    isdeg: i32,
    itdeg: i32,
    ildeg: i32,
    nsjump: usize,
    ntjump: usize,
    nljump: usize,
    ni: usize,
    userw: bool,
    rw: &mut [T],
    season: &mut [T],
    trend: &mut [T],
    work1: &mut [T],
    work2: &mut [T],
    work3: &mut [T],
    work4: &mut [T],
    work5: &mut [T],
) {
    for _ in 0..ni {
        for i in 0..n {
            work1[i] = y[i] - trend[i];
        }

        ss(
            work1, n, np, ns, isdeg, nsjump, userw, rw, work2, work3, work4, work5, season,
        );
        fts(work2, n + 2 * np, np, work3, work1);
        ess(work3, n, nl, ildeg, nljump, false, work4, work1, work5);
        for i in 0..n {
            season[i] = work2[np + i] - work1[i];
        }
        for i in 0..n {
            work1[i] = y[i] - season[i];
        }
        ess(work1, n, nt, itdeg, ntjump, userw, rw, trend, work3);
    }
}

fn rwts<T: Float>(y: &[T], n: usize, fit: &[T], rw: &mut [T]) {
    for i in 0..n {
        rw[i] = (y[i] - fit[i]).abs();
    }

    let mid1 = (n - 1) / 2;
    let mid2 = n / 2;

    rw.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

    let cmad = T::from_f64(3.0) * (rw[mid1] + rw[mid2]); // 6 * median abs resid
    let c9 = T::from_f64(0.999) * cmad;
    let c1 = T::from_f64(0.001) * cmad;

    for i in 0..n {
        let r = (y[i] - fit[i]).abs();
        if r <= c1 {
            rw[i] = T::one();
        } else if r <= c9 {
            rw[i] = pow2(T::one() - pow2(r / cmad));
        } else {
            rw[i] = T::zero();
        }
    }
}

fn ss<T: Float>(
    y: &[T],
    n: usize,
    np: usize,
    ns: usize,
    isdeg: i32,
    nsjump: usize,
    userw: bool,
    rw: &[T],
    season: &mut [T],
    work1: &mut [T],
    work2: &mut [T],
    work3: &mut [T],
    work4: &mut [T],
) {
    for j in 1..=np {
        let k = (n - j) / np + 1;

        for i in 1..=k {
            work1[i - 1] = y[(i - 1) * np + j - 1];
        }
        if userw {
            for i in 1..=k {
                work3[i - 1] = rw[(i - 1) * np + j - 1];
            }
        }
        ess(
            work1,
            k,
            ns,
            isdeg,
            nsjump,
            userw,
            work3,
            &mut work2[1..],
            work4,
        );
        let mut xs = T::zero();
        let nright = ns.min(k);
        let ok = est(
            work1,
            k,
            ns,
            isdeg,
            xs,
            &mut work2[0],
            1,
            nright,
            work4,
            userw,
            work3,
        );
        if !ok {
            work2[0] = work2[1];
        }
        xs = T::from_usize(k + 1);
        let nleft = 1.max(k as i32 - ns as i32 + 1) as usize;
        let ok = est(
            work1,
            k,
            ns,
            isdeg,
            xs,
            &mut work2[k + 1],
            nleft,
            k,
            work4,
            userw,
            work3,
        );
        if !ok {
            work2[k + 1] = work2[k];
        }
        for m in 1..=k + 2 {
            season[(m - 1) * np + j - 1] = work2[m - 1];
        }
    }
}
