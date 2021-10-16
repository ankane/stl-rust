// Ported from https://www.netlib.org/a/stl
//
// Cleveland, R. B., Cleveland, W. S., McRae, J. E., & Terpenning, I. (1990).
// STL: A Seasonal-Trend Decomposition Procedure Based on Loess.
// Journal of Official Statistics, 6(1), 3-33.

fn stl(y: &[f32], n: usize, np: usize, ns: usize, nt: usize, nl: usize, isdeg: i32, itdeg: i32, ildeg: i32, nsjump: usize, ntjump: usize, nljump: usize, ni: usize, no: usize, rw: &mut [f32], season: &mut [f32], trend: &mut [f32]) {
    let mut work1 = vec![0.0; n + 2 * np];
    let mut work2 = vec![0.0; n + 2 * np];
    let mut work3 = vec![0.0; n + 2 * np];
    let mut work4 = vec![0.0; n + 2 * np];
    let mut work5 = vec![0.0; n + 2 * np];

    let mut userw = false;
    let mut k = 0;

    // TODO add messages
    assert!(ns >= 3);
    assert!(nt >= 3);
    assert!(nl >= 3);
    assert!(np >= 2);

    assert!(isdeg == 0 || isdeg == 1);
    assert!(itdeg == 0 || itdeg == 1);
    assert!(ildeg == 0 || ildeg == 1);

    assert!(ns % 2 == 1);
    assert!(nt % 2 == 1);
    assert!(nl % 2 == 1);

    loop {
        onestp(&y, n, np, ns, nt, nl, isdeg, itdeg, ildeg, nsjump, ntjump, nljump, ni, userw, rw, season, trend, &mut work1, &mut work2, &mut work3, &mut work4, &mut work5);
        k += 1;
        if k > no {
            break;
        }
        for i in 0..n {
            work1[i] = trend[i] + season[i];
        }
        rwts(y, n, &mut work1, rw);
        userw = true;
    }

    if no <= 0 {
        for i in 0..n {
            rw[i] = 1.0;
        }
    }
}

fn ess(y: &[f32], n: usize, len: usize, ideg: i32, njump: usize, userw: bool, rw: &[f32], ys: &mut [f32], res: &mut [f32]) {
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
            let ok = est(y, n, len, ideg, i as f32, &mut ys[i - 1], nleft, nright, res, userw, rw);
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
            let ok = est(y, n, len, ideg, i as f32, &mut ys[i - 1], nleft, nright, res, userw, rw);
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
            } else if i >= n - nsh + 1 {
                nleft = n - len + 1;
                nright = n;
            } else {
                nleft = i - nsh + 1;
                nright = len + i - nsh;
            }
            let ok = est(y, n, len, ideg, i as f32, &mut ys[i - 1], nleft, nright, res, userw, rw);
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
            let ok = est(y, n, len, ideg, n as f32, &mut ys[n - 1], nleft, nright, res, userw, rw);
            if !ok {
                ys[n - 1] = y[n - 1];
                if k != n - 1 {
                    let delta = (ys[n - 1] - ys[k - 1]) / ((n - k) as f32);
                    for j in k + 1..=n - 1 {
                        ys[j - 1] = ys[k - 1] + delta * ((j - k) as f32);
                    }
                }
            }
        }
    }
}

fn est(y: &[f32], n: usize, len: usize, ideg: i32, xs: f32, ys: &mut f32, nleft: usize, nright: usize, w: &mut [f32], userw: bool, rw: &[f32]) -> bool {
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
                w[j - 1] = (1.0 - (r / h).powi(3)).powi(3);
            }
            if userw {
                w[j - 1] *= rw[j - 1];
            }
            a += w[j - 1];
        }
    }

    if a <= 0.0 {
        false
    } else { // weighted least squares
        for j in nleft..=nright { // make sum of w(j) == 1
            w[j - 1] /= a;
        }

        if h > 0.0 && ideg > 0 { // use linear fit
            let mut a = 0.0;
            for j in nleft..=nright { // weighted center of x values
                a += w[j - 1] * (j as f32);
            }
            let mut b = xs - a;
            let mut c = 0.0;
            for j in nleft..=nright {
                c += w[j - 1] * ((j as f32) - a).powi(2);
            }
            if c.sqrt() > 0.001 * range {
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
    let mut v = 0.0;

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

fn onestp(y: &[f32], n: usize, np: usize, ns: usize, nt: usize, nl: usize, isdeg: i32, itdeg: i32, ildeg: i32, nsjump: usize, ntjump: usize, nljump: usize, ni: usize, userw: bool, rw: &mut [f32], season: &mut [f32], trend: &mut [f32], work1: &mut [f32], work2: &mut [f32], work3: &mut [f32], work4: &mut [f32], work5: &mut [f32]) {
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

fn rwts(y: &[f32], n: usize, fit: &[f32], rw: &mut [f32]) {
    for i in 0..n {
        rw[i] = (y[i] - fit[i]).abs();
    }

    let mid1 = (n - 1) / 2;
    let mid2 = n / 2;

    // sort
    rw.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let cmad = 3.0 * (rw[mid1] + rw[mid2]); // 6 * median abs resid
    let c9 = 0.999 * cmad;
    let c1 = 0.001 * cmad;

    for i in 0..n {
        let r = (y[i] - fit[i]).abs();
        if r <= c1 {
            rw[i] = 1.0;
        } else if r <= c9 {
            rw[i] = (1.0 - (r / cmad).powi(2)).powi(2);
        } else {
            rw[i] = 0.0;
        }
    }
}

fn ss(y: &[f32], n: usize, np: usize, ns: usize, isdeg: i32, nsjump: usize, userw: bool, rw: &[f32], season: &mut [f32], work1: &mut [f32], work2: &mut [f32], work3: &mut [f32], work4: &mut [f32]) {
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
        let mut xs = 0.0;
        let nright = ns.min(k);
        let ok = est(work1, k, ns, isdeg, xs, &mut work2[0], 1, nright, work4, userw, work3);
        if !ok {
            work2[0] = work2[1];
        }
        xs = (k + 1) as f32;
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

#[derive(Debug)]
pub struct StlParams {
    ns: Option<usize>,
    nt: Option<usize>,
    nl: Option<usize>,
    isdeg: i32,
    itdeg: i32,
    ildeg: Option<i32>,
    nsjump: Option<usize>,
    ntjump: Option<usize>,
    nljump: Option<usize>,
    ni: Option<usize>,
    no: Option<usize>,
    robust: bool
}

#[derive(Debug)]
pub struct StlResult {
    seasonal: Vec<f32>,
    trend: Vec<f32>,
    remainder: Vec<f32>,
    weights: Vec<f32>
}

pub fn params() -> StlParams {
    StlParams {
        ns: None,
        nt: None,
        nl: None,
        isdeg: 0,
        itdeg: 1,
        ildeg: None,
        nsjump: None,
        ntjump: None,
        nljump: None,
        ni: None,
        no: None,
        robust: false
    }
}

impl StlParams {
    pub fn seasonal_length(&mut self, ns: usize) -> &mut Self {
        self.ns = Some(ns);
        self
    }

    pub fn trend_length(&mut self, nt: usize) -> &mut Self {
        self.nt = Some(nt);
        self
    }

    pub fn low_pass_length(&mut self, nl: usize) -> &mut Self {
        self.nl = Some(nl);
        self
    }

    pub fn seasonal_degree(&mut self, isdeg: i32) -> &mut Self {
        self.isdeg = isdeg;
        self
    }

    pub fn trend_degree(&mut self, itdeg: i32) -> &mut Self {
        self.itdeg = itdeg;
        self
    }

    pub fn low_pass_degree(&mut self, ildeg: i32) -> &mut Self {
        self.ildeg = Some(ildeg);
        self
    }

    pub fn seasonal_jump(&mut self, nsjump: usize) -> &mut Self {
        self.nsjump = Some(nsjump);
        self
    }

    pub fn trend_jump(&mut self, ntjump: usize) -> &mut Self {
        self.ntjump = Some(ntjump);
        self
    }

    pub fn low_pass_jump(&mut self, nljump: usize) -> &mut Self {
        self.nljump = Some(nljump);
        self
    }

    pub fn inner_loops(&mut self, ni: usize) -> &mut Self {
        self.ni = Some(ni);
        self
    }

    pub fn outer_loops(&mut self, no: usize) -> &mut Self {
        self.no = Some(no);
        self
    }

    pub fn robust(&mut self, robust: bool) -> &mut Self {
        self.robust = robust;
        self
    }

    pub fn fit(&self, y: &[f32], np: usize) -> StlResult {
        let n = y.len();

        assert!(n >= np * 2, "series has less than two periods");

        let ns = self.ns.unwrap_or(np);

        let isdeg = self.isdeg;
        let itdeg = self.itdeg;

        let mut rw = vec![0.0; n];
        let mut season = vec![0.0; n];
        let mut trend = vec![0.0; n];

        let ildeg = self.ildeg.unwrap_or(itdeg);
        let mut newns = ns.max(3);
        if newns % 2 == 0 {
            newns += 1;
        }

        let newnp = np.max(2);
        let mut nt = ((1.5 * newnp as f32) / (1.0 - 1.5 / newns as f32)).ceil() as usize;
        nt = self.nt.unwrap_or(nt);
        nt = nt.max(3);
        if nt % 2 == 0 {
            nt += 1;
        }

        let mut nl = self.nl.unwrap_or(newnp);
        if nl % 2 == 0 && self.nl.is_none() {
            nl += 1;
        }

        let ni = self.ni.unwrap_or(if self.robust { 1 } else { 2 });
        let no = self.no.unwrap_or(if self.robust { 15 } else { 0 });

        let nsjump = self.nsjump.unwrap_or(((newns as f32) / 10.0).ceil() as usize);
        let ntjump = self.ntjump.unwrap_or(((nt as f32) / 10.0).ceil() as usize);
        let nljump = self.nljump.unwrap_or(((nl as f32) / 10.0).ceil() as usize);

        stl(y, n, newnp, newns, nt, nl, isdeg, itdeg, ildeg, nsjump, ntjump, nljump, ni, no, &mut rw, &mut season, &mut trend);

        let mut remainder = Vec::with_capacity(n);
        for i in 0..n {
            remainder.push(y[i] - season[i] - trend[i]);
        }

        StlResult {
            seasonal: season,
            trend,
            remainder,
            weights: rw
        }
    }
}

impl StlResult {
    pub fn seasonal(&self) -> &Vec<f32> {
        &self.seasonal
    }

    pub fn trend(&self) -> &Vec<f32> {
        &self.trend
    }

    pub fn remainder(&self) -> &Vec<f32> {
        &self.remainder
    }

    pub fn weights(&self) -> &Vec<f32> {
        &self.weights
    }
}

#[cfg(test)]
mod tests {
    fn assert_elements_in_delta(exp: &[f32], act: &[f32]) {
        assert_eq!(exp.len(), act.len());
        for i in 0..exp.len() {
            assert!((exp[i] - act[i]).abs() < 0.001);
        }
    }

    fn generate_series() -> Vec<f32> {
        return vec![
            5.0, 9.0, 2.0, 9.0, 0.0, 6.0, 3.0, 8.0, 5.0, 8.0,
            7.0, 8.0, 8.0, 0.0, 2.0, 5.0, 0.0, 5.0, 6.0, 7.0,
            3.0, 6.0, 1.0, 4.0, 4.0, 4.0, 3.0, 7.0, 5.0, 8.0
        ];
    }

    #[test]
    fn test_works() {
        let series = generate_series();
        let result = crate::params().fit(&series, 7);
        assert_elements_in_delta(&[0.36926576, 0.75655484, -1.3324139, 1.9553658, -0.6044802], &result.seasonal()[..5]);
        assert_elements_in_delta(&[4.804099, 4.9097075, 5.015316, 5.16045, 5.305584], &result.trend()[..5]);
        assert_elements_in_delta(&[-0.17336464, 3.3337379, -1.6829021, 1.8841844, -4.7011037], &result.remainder()[..5]);
        assert_elements_in_delta(&[1.0, 1.0, 1.0, 1.0, 1.0], &result.weights()[..5]);
    }

    #[test]
    fn test_robust() {
        let series = generate_series();
        let result = crate::params().robust(true).fit(&series, 7);
        assert_elements_in_delta(&[0.14922355, 0.47939026, -1.833231, 1.7411387, 0.8200711], &result.seasonal()[..5]);
        assert_elements_in_delta(&[5.397365, 5.4745436, 5.5517216, 5.6499176, 5.748114], &result.trend()[..5]);
        assert_elements_in_delta(&[-0.5465884, 3.0460663, -1.7184906, 1.6089439, -6.5681853], &result.remainder()[..5]);
        assert_elements_in_delta(&[0.99374926, 0.8129377, 0.9385952, 0.9458036, 0.29742217], &result.weights()[..5]);
    }

    #[test]
    #[should_panic(expected = "series has less than two periods")]
    fn test_too_few_periods() {
        let series = generate_series();
        crate::params().fit(&series, 16);
    }
}
