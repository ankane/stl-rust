// Ported from https://www.netlib.org/a/stl
//
// Cleveland, R. B., Cleveland, W. S., McRae, J. E., & Terpenning, I. (1990).
// STL: A Seasonal-Trend Decomposition Procedure Based on Loess.
// Journal of Official Statistics, 6(1), 3-33.

#![allow(clippy::too_many_arguments)]

#[cfg(feature = "std")]
fn sqrt(x: f32) -> f32 {
    x.sqrt()
}

#[cfg(not(feature = "std"))]
use core::f32::math::sqrt;

fn pow2(x: f32) -> f32 {
    x * x
}

fn pow3(x: f32) -> f32 {
    x * x * x
}

pub fn stl(
    y: &[f32],
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
    rw: &mut [f32],
    season: &mut [f32],
    trend: &mut [f32],
    work: &mut [f32],
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
            *v = 1.0;
        }
    }
}

fn ess(
    y: &[f32],
    n: usize,
    len: usize,
    ideg: i32,
    njump: usize,
    userw: bool,
    rw: &[f32],
    ys: &mut [f32],
    res: &mut [f32],
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
                i as f32,
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
                i as f32,
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
                i as f32,
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
            let delta = (ys[i + newnj - 1] - ys[i - 1]) / (newnj as f32);
            for j in i + 1..=i + newnj - 1 {
                ys[j - 1] = ys[i - 1] + delta * ((j - i) as f32);
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
                n as f32,
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
                let delta = (ys[n - 1] - ys[k - 1]) / ((n - k) as f32);
                for j in k + 1..=n - 1 {
                    ys[j - 1] = ys[k - 1] + delta * ((j - k) as f32);
                }
            }
        }
    }
}

fn est(
    y: &[f32],
    n: usize,
    len: usize,
    ideg: i32,
    xs: f32,
    ys: &mut f32,
    nleft: usize,
    nright: usize,
    w: &mut [f32],
    userw: bool,
    rw: &[f32],
) -> bool {
    let range = (n as f32) - 1.0;
    let mut h = (xs - (nleft as f32)).max((nright as f32) - xs);

    if len > n {
        h += ((len - n) / 2) as f32;
    }

    let h9 = 0.999 * h;
    let h1 = 0.001 * h;

    // compute weights
    let mut a = 0.0;
    for j in nleft..=nright {
        w[j - 1] = 0.0;
        let r = ((j as f32) - xs).abs();
        if r <= h9 {
            if r <= h1 {
                w[j - 1] = 1.0;
            } else {
                w[j - 1] = pow3(1.0 - pow3(r / h));
            }
            if userw {
                w[j - 1] *= rw[j - 1];
            }
            a += w[j - 1];
        }
    }

    if a <= 0.0 {
        false
    } else {
        // weighted least squares
        for j in nleft..=nright {
            // make sum of w(j) == 1
            w[j - 1] /= a;
        }

        if h > 0.0 && ideg > 0 {
            // use linear fit
            let mut a = 0.0;
            for j in nleft..=nright {
                // weighted center of x values
                a += w[j - 1] * (j as f32);
            }
            let mut b = xs - a;
            let mut c = 0.0;
            for j in nleft..=nright {
                c += w[j - 1] * pow2((j as f32) - a);
            }
            if sqrt(c) > 0.001 * range {
                b /= c;

                // points are spread out enough to compute slope
                for j in nleft..=nright {
                    w[j - 1] *= b * ((j as f32) - a) + 1.0;
                }
            }
        }

        *ys = 0.0;
        for j in nleft..=nright {
            *ys += w[j - 1] * y[j - 1];
        }

        true
    }
}

fn fts(x: &[f32], n: usize, np: usize, trend: &mut [f32], work: &mut [f32]) {
    ma(x, n, np, trend);
    ma(trend, n - np + 1, np, work);
    ma(work, n - 2 * np + 2, 3, trend);
}

fn ma(x: &[f32], n: usize, len: usize, ave: &mut [f32]) {
    let newn = n - len + 1;
    let flen = len as f32;

    // get the first average
    let mut v: f32 = x.iter().take(len).sum();
    ave[0] = v / flen;

    if newn > 1 {
        for (m, aj) in ave.iter_mut().take(newn).skip(1).enumerate() {
            // window down the array
            v = v - x[m] + x[len + m];
            *aj = v / flen;
        }
    }
}

fn onestp(
    y: &[f32],
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
    rw: &mut [f32],
    season: &mut [f32],
    trend: &mut [f32],
    work1: &mut [f32],
    work2: &mut [f32],
    work3: &mut [f32],
    work4: &mut [f32],
    work5: &mut [f32],
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

fn rwts(y: &[f32], n: usize, fit: &[f32], rw: &mut [f32]) {
    for i in 0..n {
        rw[i] = (y[i] - fit[i]).abs();
    }

    let mid1 = (n - 1) / 2;
    let mid2 = n / 2;

    rw.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

    let cmad = 3.0 * (rw[mid1] + rw[mid2]); // 6 * median abs resid
    let c9 = 0.999 * cmad;
    let c1 = 0.001 * cmad;

    for i in 0..n {
        let r = (y[i] - fit[i]).abs();
        if r <= c1 {
            rw[i] = 1.0;
        } else if r <= c9 {
            rw[i] = pow2(1.0 - pow2(r / cmad));
        } else {
            rw[i] = 0.0;
        }
    }
}

fn ss(
    y: &[f32],
    n: usize,
    np: usize,
    ns: usize,
    isdeg: i32,
    nsjump: usize,
    userw: bool,
    rw: &[f32],
    season: &mut [f32],
    work1: &mut [f32],
    work2: &mut [f32],
    work3: &mut [f32],
    work4: &mut [f32],
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
        let mut xs = 0.0;
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
        xs = (k + 1) as f32;
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
