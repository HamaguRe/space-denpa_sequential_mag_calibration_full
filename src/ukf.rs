//! UKF版の逐次推定アルゴリズム
//! 
//! 参考文献
//! 1. John L. Crassdis, Kok-Lam Lai, Richard R. Harman, 
//!    "Real-Time Attitude-Independent Three-Axis Magnetometer Calibration", 
//!    Journal of Guidance, Control, and Dynamics, Vol.28, No.1, 2005.

use super::{SVector, SMatrix, NOISE_VAR, FIELD_NORM};

type Vector3 = SVector<f64, 3>;
type Vector9 = SVector<f64, 9>;
type Matrix3x3 = SMatrix<f64, 3, 3>;
type Matrix9x9 = SMatrix<f64, 9, 9>;

/// 状態変数の個数
const N_X: usize = 9;

/// シグマポイントの拡げ具合を調整するパラメータ
const ALPHA: f64 = 0.1;

/// 分布の事前情報を組み込むためのパラメータ．正規分布なら2が最適[Ref.1]．
const BETA: f64 = 2.0;  // 大きくすると較正した値の描く球が大きくなる

const KAPPA: f64 = 3.0 - N_X as f64;
const LAMBDA: f64 = ALPHA * ALPHA * (N_X as f64 + KAPPA) - N_X as f64;


pub struct UKF {
    pub xhat: Vector9,  // 状態変数（先頭3要素がバイアス，その後ろ3要素がスケールファクター，最後の3要素は非直交性補正係数）
    p: Matrix9x9,  // 誤差共分散行列
}

impl UKF {
    pub fn new() -> Self {
        let mut p_init = Matrix9x9::zeros();
        for i in 0..3 {
            p_init[(i, i)]     = 0.1;    // バイアス
            p_init[(i+3, i+3)] = 0.001;  // スケールファクター
            p_init[(i+6, i+6)] = 0.0001; // 非直交性補正係数
        }

        Self {
            xhat: Vector9::zeros(),
            p: p_init,
        }
    }

    pub fn update(&mut self, mag: Vector3) {
        // --- シグマポイントを計算 ---
        let l = self.p.cholesky().unwrap().unpack();
        let mut sigma_points = [Vector9::zeros(); 2*N_X + 1];
        sigma_points[0] = self.xhat;
        let coef = (N_X as f64 + LAMBDA).sqrt();
        for i in 0..N_X {
            let tmp = coef * l.column(i);
            sigma_points[i + 1]       = self.xhat + tmp;
            sigma_points[i + 1 + N_X] = self.xhat - tmp;
        }

        // 出力のシグマポイントの更新
        let mut z_sigma_points = [0.0; 2*N_X + 1];
        for i in 0..(2*N_X + 1) {
            z_sigma_points[i] = h(sigma_points[i], mag);
        }

        // 事前出力推定値
        let mut zhat = 0.0;
        for i in 0..(2*N_X + 1) {
            zhat += w_mean(i) * z_sigma_points[i];
        }
        
        // 事前出力誤差共分散行列
        let mut p_zz = 0.0;
        for i in 0..(2*N_X + 1) {
            let tmp = z_sigma_points[i] - zhat;
            p_zz += w_cov(i) * tmp * tmp;
        }

        // 事前状態・出力誤差共分散行列
        let mut p_xz = Vector9::zeros();
        for i in 0..(2*N_X + 1) {
            p_xz += w_cov(i) * (sigma_points[i] - self.xhat) * (z_sigma_points[i] - zhat);
        }

        
        // ノイズの共分散行列(Ref1 - eq.(5c))は定数になるので計算後の値を埋め込んである
        let tmp = self.calibrate(mag);
        let sigma_square = (4.0 * tmp.transpose() * NOISE_VAR * tmp)[0] + 2.0 * (3.0 * NOISE_VAR * NOISE_VAR);

        // --- フィルタリングステップ ---
        // カルマンゲイン
        let gain = p_xz / (p_zz + sigma_square);
        // 状態推定値
        let z = mag.norm_squared() - FIELD_NORM * FIELD_NORM;
        self.xhat = self.xhat + gain * (z - zhat);
        // 事後誤差共分散行列
        self.p = self.p - gain * (p_zz + sigma_square) * gain.transpose();
        // ----------------------------------- //
    }

    /// 推定したバイアスとスケールファクタを用いて較正した地磁気計測値を返す
    /// 
    /// * mag: 生の地磁気計測値 [x, y, z]
    pub fn calibrate(&mut self, mag: Vector3) -> Vector3 {
        let (b, d) = get_bias_scale_factor(self.xhat);

        (mag - b) + d * mag
    }
}

/// 状態変数xからバイアスとスケールファクタを取り出す．
fn get_bias_scale_factor(x: Vector9) -> (Vector3, Matrix3x3) {
    let mut d = Matrix3x3::zeros();  // 対称行列
    for i in 0..3 {
        d[(i, i)] = x[i + 3]
    }
    d[(0, 1)] = x[6];
    d[(0, 2)] = x[7];
    d[(1, 2)] = x[8];
    d[(1, 0)] = x[6];
    d[(2, 0)] = x[7];
    d[(2, 1)] = x[8];

    (x.fixed_rows::<3>(0).into_owned(), d)
}

fn h(x: Vector9, mag: Vector3) -> f64 {
    let (b, d) = get_bias_scale_factor(x);

    (-mag.transpose() * (2.0 * d + d * d) * mag + 2.0 * mag.transpose() * (Matrix3x3::identity() + d) * b)[0] - b.norm_squared()
}

// なんか微妙な感じするけど，2n+1個の配列を持つよりはこっちのほうが良い
/// 重み係数w_cov(i)を計算
fn w_cov(i: usize) -> f64 {
    if i == 0 {
        (LAMBDA / (N_X as f64 + LAMBDA)) + (1.0 - ALPHA * ALPHA + BETA)
    } else {
        1.0 / (2.0 * (N_X as f64 + LAMBDA))
    }
}

/// 重み係数w_mean(i)を計算
fn w_mean(i: usize) -> f64 {
    if i == 0 {
        LAMBDA / (N_X as f64 + LAMBDA)
    } else {
        1.0 / (2.0 * (N_X as f64 + LAMBDA))
    }
}
