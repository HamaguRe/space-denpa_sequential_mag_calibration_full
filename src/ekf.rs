//! EKF版の逐次推定アルゴリズム

use super::{SVector, SMatrix, NOISE_VAR, FIELD_NORM};

type Scalar = SVector<f64, 1>;  // 1x1行列とf64を足せないからスカラー型を作る
type Vector3 = SVector<f64, 3>;
type Vector9 = SVector<f64, 9>;
type Matrix1x6 = SMatrix<f64, 1, 6>;
type Matrix1x9 = SMatrix<f64, 1, 9>;
type Matrix3x3 = SMatrix<f64, 3, 3>;
type Matrix6x6 = SMatrix<f64, 6, 6>;
type Matrix9x9 = SMatrix<f64, 9, 9>;

pub struct EKF {
    pub xhat: Vector9,  // 状態変数（先頭3要素がバイアス，その後ろ3要素がスケールファクター，最後の3要素は非直交性補正係数）
    p: Matrix9x9, // 誤差共分散行列
}

impl EKF {
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
        #[allow(non_snake_case)]
        let H = self.jacobian_h(mag);

        // ノイズの共分散行列(Ref1 - eq.(5c))は定数になるので計算後の値を埋め込んである
        let tmp = self.calibrate(mag);
        let sigma_square = 4.0 * tmp.transpose() * NOISE_VAR * tmp + Scalar::new(2.0 * (3.0 * NOISE_VAR * NOISE_VAR));

        let coef: f64 = 1.0 / (H * self.p * H.transpose() + sigma_square)[0];  // スカラーになるので逆行列計算は不要
        let gain = self.p * H.transpose() * coef;
        let z =  Scalar::new(mag.norm_squared() - FIELD_NORM * FIELD_NORM);
        self.xhat = self.xhat + gain * (z - self.h(mag));
        self.p = (Matrix9x9::identity() - gain * H) * self.p;
    }

    fn get_bias_scale_factor(&self) -> (Vector3, Matrix3x3) {
        let mut d = Matrix3x3::zeros();
        for i in 0..3 {
            d[(i, i)] = self.xhat[i + 3]
        }
        d[(0, 1)] = self.xhat[6];
        d[(0, 2)] = self.xhat[7];
        d[(1, 2)] = self.xhat[8];
        d[(1, 0)] = self.xhat[6];
        d[(2, 0)] = self.xhat[7];
        d[(2, 1)] = self.xhat[8];

        (self.xhat.fixed_rows::<3>(0).into_owned(), d)
    }

    /// 推定したバイアスとスケールファクタを用いて較正した地磁気計測値を返す
    /// 
    /// * mag: 生の地磁気計測値 [x, y, z]
    pub fn calibrate(&self, mag: Vector3) -> Vector3 {
        let (b, d) = self.get_bias_scale_factor();

        (mag - b) + d * mag
    }

    fn h(&self, mag: Vector3) -> Scalar {
        let (b, d) = self.get_bias_scale_factor();
        -mag.transpose() * (2.0 * d + d * d) * mag + 2.0 * mag.transpose() * (Matrix3x3::identity() + d) * b - Scalar::new(b.norm_squared())
    }

    fn jacobian_h(&self, mag: Vector3) -> Matrix1x9 {
        let (b, d) = self.get_bias_scale_factor();

        #[allow(non_snake_case)]
        let B1 = mag[0];
        #[allow(non_snake_case)]
        let B2 = mag[1];
        #[allow(non_snake_case)]
        let B3 = mag[2];
        let s = Matrix1x6::new(B1*B1, B2*B2, B3*B3, 2.0*B1*B2, 2.0*B1*B3, 2.0*B2*B3);
        let mut de_dd = Matrix6x6::zeros();
        // 対角成分
        de_dd[(0, 0)] = 2.0 * (1.0 + d[(0, 0)]);
        de_dd[(1, 1)] = 2.0 * (1.0 + d[(1, 1)]);
        de_dd[(2, 2)] = 2.0 * (1.0 + d[(2, 2)]);
        de_dd[(3, 3)] = 2.0 + d[(0, 0)] + d[(1, 1)];
        de_dd[(4, 4)] = 2.0 + d[(0, 0)] + d[(2, 2)];
        de_dd[(5, 5)] = 2.0 + d[(1, 1)] + d[(2, 2)];
        // 上三角成分
        de_dd[(0, 3)] = 2.0 * d[(0, 1)];
        de_dd[(0, 4)] = 2.0 * d[(0, 2)];
        de_dd[(1, 3)] = 2.0 * d[(0, 1)];
        de_dd[(1, 5)] = 2.0 * d[(1, 2)];
        de_dd[(2, 4)] = 2.0 * d[(0, 2)];
        de_dd[(2, 5)] = 2.0 * d[(1, 2)];
        de_dd[(3, 4)] = d[(1, 2)];
        de_dd[(3, 5)] = d[(0, 2)];
        de_dd[(4, 5)] = d[(0, 1)];
        // 下三角成分
        de_dd[(3, 0)] = d[(0, 1)];
        de_dd[(3, 1)] = d[(0, 1)];
        de_dd[(4, 0)] = d[(0, 2)];
        de_dd[(4, 2)] = d[(0, 2)];
        de_dd[(4, 3)] = d[(1, 2)];
        de_dd[(5, 1)] = d[(1, 2)];
        de_dd[(5, 2)] = d[(1, 2)];
        de_dd[(5, 3)] = d[(0, 2)];
        de_dd[(5, 4)] = d[(0, 1)];

        let b1 = self.xhat[0];
        let b2 = self.xhat[1];
        let b3 = self.xhat[2];
        let j = Matrix1x6::new(
            B1*b1, B2*b2, B3*b3,
            B1*b2 + B2*b1, B1*b3 + B3*b1, B2*b3 + B3*b2,
        );
        
        let tmp1x3 = 2.0 * mag.transpose() * (Matrix3x3::identity() + d) -  2.0 * b.transpose();
        let tmp1x6 = -s * de_dd + 2.0 * j;
        let mut dh_dx = Matrix1x9::zeros();
        // 行列を結合する方法無いんだっけ
        for i in 0..3 {
            dh_dx[(0, i)] = tmp1x3[(0, i)];
        }
        for i in 0..6 {
            dh_dx[(0, i + 3)] = tmp1x6[(0, i)];
        }
        dh_dx
    }
}