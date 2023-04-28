// Ported from https://www.netlib.org/a/stl
//
// Cleveland, R. B., Cleveland, W. S., McRae, J. E., & Terpenning, I. (1990).
// STL: A Seasonal-Trend Decomposition Procedure Based on Loess.
// Journal of Official Statistics, 6(1), 3-33.

use crate::Float;

pub fn stl<F: Float>(y: &[F], n: usize, np: usize, ns: usize, nt: usize, nl: usize, isdeg: i32, itdeg: i32, ildeg: i32, nsjump: usize, ntjump: usize, nljump: usize, ni: usize, no: usize, rw: &mut [F], season: &mut [F], trend: &mut [F]) {
    let mut work1: Vec<F> = vec![F::zero(); n + 2 * np];
    let mut work2: Vec<F> = vec![F::zero(); n + 2 * np];
    let mut work3: Vec<F> = vec![F::zero(); n + 2 * np];
    let mut work4: Vec<F> = vec![F::zero(); n + 2 * np];
    let mut work5: Vec<F> = vec![F::zero(); n + 2 * np];

    let mut userw = false;
    let mut k = 0;

    loop {
        onestp(y, n, np, ns, nt, nl, isdeg, itdeg, ildeg, nsjump, ntjump, nljump, ni, userw, rw, season, trend, &mut work1, &mut work2, &mut work3, &mut work4, &mut work5);
        k += 1;
        if k > no {
            break;
        }
        for i in 0..n {
            work1[i] = trend[i] + season[i];
        }
        rwts(y, n, &work1, rw);
        userw = true;
    }

    if no == 0 {
        for i in 0..n {
            rw[i] = F::one();
        }
    }
}

fn ess<F: Float>(y: &[F], n: usize, len: usize, ideg: i32, njump: usize, userw: bool, rw: &[F], ys: &mut [F], res: &mut [F]) {
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
            let ok = est(y, n, len, ideg, F::from_usize(i).unwrap(), &mut ys[i - 1], nleft, nright, res, userw, rw);
            if !ok {
                ys[i - 1] = y[i - 1];
            }
            i += newnj;
        }
    } else if newnj == 1 { // newnj equal to one, len less than n
        let nsh = (len + 1) / 2;
        nleft = 1;
        nright = len;
        for i in 1..=n { // fitted value at i
            if i > nsh && nright != n {
                nleft += 1;
                nright += 1;
            }
            let ok = est(y, n, len, ideg, F::from_usize(i).unwrap(), &mut ys[i - 1], nleft, nright, res, userw, rw);
            if !ok {
                ys[i - 1] = y[i - 1];
            }
        }
    } else { // newnj greater than one, len less than n
        let nsh = (len + 1) / 2;
        let mut i = 1;
        while i <= n { // fitted value at i
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
            let ok = est(y, n, len, ideg, F::from_usize(i).unwrap(), &mut ys[i - 1], nleft, nright, res, userw, rw);
            if !ok {
                ys[i - 1] = y[i - 1];
            }
            i += newnj;
        }
    }

    if newnj != 1 {
        let mut i = 1;
        while i <= n - newnj {
            let delta = (ys[i + newnj - 1] - ys[i - 1]) / F::from_usize(newnj).unwrap();
            for j in i + 1..=i + newnj - 1 {
                ys[j - 1] = ys[i - 1] + delta * F::from_usize(j - i).unwrap();
            }
            i += newnj;
        }
        let k = ((n - 1) / newnj) * newnj + 1;
        if k != n {
            let ok = est(y, n, len, ideg, F::from_usize(n).unwrap(), &mut ys[n - 1], nleft, nright, res, userw, rw);
            if !ok {
                ys[n - 1] = y[n - 1];
                if k != n - 1 {
                    let delta = (ys[n - 1] - ys[k - 1]) / F::from_usize(n - k).unwrap();
                    for j in k + 1..=n - 1 {
                        ys[j - 1] = ys[k - 1] + delta * F::from_usize(j - k).unwrap();
                    }
                }
            }
        }
    }
}

fn est<F: Float>(y: &[F], n: usize, len: usize, ideg: i32, xs: F, ys: &mut F, nleft: usize, nright: usize, w: &mut [F], userw: bool, rw: &[F]) -> bool {
    let range = F::from_usize(n).unwrap() - F::one();
    let mut h = (xs - F::from_usize(nleft).unwrap()).max( F::from_usize(nright).unwrap() - xs);

    if len > n {
        h += F::from_usize((len - n) / 2).unwrap();
    }

    let h9 = <F as From<_>>::from(0.999) * h;
    let h1 = <F as From<_>>::from(0.001) * h;

    // compute weights
    let mut a = F::zero();
    for j in nleft..=nright {
        w[j - 1] = F::zero();
        let r = (F::from_usize(j).unwrap() - xs).abs();
        if r <= h9 {
            if r <= h1 {
                w[j - 1] = F::one();
            } else {
                w[j - 1] = (F::one() - (r / h).powi(3)).powi(3);
            }
            if userw {
                w[j - 1] *= rw[j - 1];
            }
            a += w[j - 1];
        }
    }

    if a <= F::zero() {
        false
    } else { // weighted least squares
        for j in nleft..=nright { // make sum of w(j) == 1
            w[j - 1] /= a;
        }

        if h > F::zero() && ideg > 0 { // use linear fit
            let mut a = F::zero();
            for j in nleft..=nright { // weighted center of x values
                a += w[j - 1] * F::from_usize(j).unwrap();
            }
            let mut b = xs - a;
            let mut c = F::zero();
            for j in nleft..=nright {
                c += w[j - 1] * (F::from_usize(j).unwrap() - a).powi(2);
            }
            if c.sqrt() > <F as From<_>>::from(0.001) * range {
                b /= c;

                // points are spread out enough to compute slope
                for j in nleft..=nright {
                    w[j - 1] *= b * (F::from_usize(j).unwrap() - a) + F::one();
                }
            }
        }

        *ys = F::zero();
        for j in nleft..=nright {
            *ys += w[j - 1] * y[j - 1];
        }

        true
    }
}

fn fts<F: Float>(x: &[F], n: usize, np: usize, trend: &mut [F], work: &mut [F]) {
    ma(x, n, np, trend);
    ma(trend, n - np + 1, np, work);
    ma(work, n - 2 * np + 2, 3, trend);
}

fn ma<F: Float>(x: &[F], n: usize, len: usize, ave: &mut [F]) {
    let newn = n - len + 1;
    let flen = F::from_usize(len).unwrap();
    let mut v = F::zero();

    // get the first average
    for i in 0..len {
        v += x[i];
    }

    ave[0] = v / flen;
    if newn > 1 {
        let mut k = len;
        let mut m = 0;
        for j in 1..newn {
            // window down the array
            v = v - x[m] + x[k];
            ave[j] = v / flen;
            k += 1;
            m += 1;
        }
    }
}

fn onestp<F: Float>(y: &[F], n: usize, np: usize, ns: usize, nt: usize, nl: usize, isdeg: i32, itdeg: i32, ildeg: i32, nsjump: usize, ntjump: usize, nljump: usize, ni: usize, userw: bool, rw: &mut [F], season: &mut [F], trend: &mut [F], work1: &mut [F], work2: &mut [F], work3: &mut [F], work4: &mut [F], work5: &mut [F]) {
    for _ in 0..ni {
        for i in 0..n {
            work1[i] = y[i] - trend[i];
        }

        ss(work1, n, np, ns, isdeg, nsjump, userw, rw, work2, work3, work4, work5, season);
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

fn rwts<F: Float>(y: &[F], n: usize, fit: &[F], rw: &mut [F]) {
    for i in 0..n {
        rw[i] = (y[i] - fit[i]).abs();
    }

    let mid1 = (n - 1) / 2;
    let mid2 = n / 2;

    rw.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

    let cmad = <F as From<_>>::from(3.0) * (rw[mid1] + rw[mid2]); // 6 * median abs resid
    let c9 = <F as From<_>>::from(0.999) * cmad;
    let c1 = <F as From<_>>::from(0.001) * cmad;

    for i in 0..n {
        let r = (y[i] - fit[i]).abs();
        if r <= c1 {
            rw[i] = F::one();
        } else if r <= c9 {
            rw[i] = (F::one() - (r / cmad).powi(2)).powi(2);
        } else {
            rw[i] = F::zero();
        }
    }
}

fn ss<F: Float>(y: &[F], n: usize, np: usize, ns: usize, isdeg: i32, nsjump: usize, userw: bool, rw: &[F], season: &mut [F], work1: &mut [F], work2: &mut [F], work3: &mut [F], work4: &mut [F]) {
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
        ess(work1, k, ns, isdeg, nsjump, userw, work3, &mut work2[1..], work4);
        let mut xs: F = F::zero();
        let nright = ns.min(k);
        let ok = est(work1, k, ns, isdeg, xs, &mut work2[0], 1, nright, work4, userw, work3);
        if !ok {
            work2[0] = work2[1];
        }
        xs = F::from_usize(k + 1).unwrap();
        let nleft = 1.max(k as i32 - ns as i32 + 1) as usize;
        let ok = est(work1, k, ns, isdeg, xs, &mut work2[k + 1], nleft, k, work4, userw, work3);
        if !ok {
            work2[k + 1] = work2[k];
        }
        for m in 1..=k + 2 {
            season[(m - 1) * np + j - 1] = work2[m - 1];
        }
    }
}
