use num_traits::{FromPrimitive, MulAdd, Num, ToPrimitive};
use crate::sos::SisoSosFilter;
mod butter_sos_tables_2;
mod butter_sos_tables_4;
mod butter_sos_tables_6;

pub fn butter2<T>(cutoff_ratio: f64) -> Result<SisoSosFilter<1, T>, &'static str>
where
    T: Num + Copy + MulAdd<Output = T> + FromPrimitive + ToPrimitive,
{
    // Convert SOS tables from static arrays to slices because new_interpolated expects slices.
    let sos_tables = butter_sos_tables_2::SOS_TABLES
        .each_ref()
        .map(|sec| sec.each_ref().map(|coeffs| &coeffs[..]));
    SisoSosFilter::new_interpolated(
        cutoff_ratio,
        &butter_sos_tables_2::LOG10_CUTOFF_RATIOS,
        sos_tables,
    )
}

pub fn butter4<T>(cutoff_ratio: f64) -> Result<SisoSosFilter<2, T>, &'static str>
where
    T: Num + Copy + MulAdd<Output = T> + FromPrimitive + ToPrimitive,
{
    // Convert SOS tables from static arrays to slices because new_interpolated expects slices.
    let sos_tables = butter_sos_tables_4::SOS_TABLES
        .each_ref()
        .map(|sec| sec.each_ref().map(|coeffs| &coeffs[..]));
    SisoSosFilter::new_interpolated(
        cutoff_ratio,
        &butter_sos_tables_4::LOG10_CUTOFF_RATIOS,
        sos_tables,
    )
}

pub fn butter6<T>(cutoff_ratio: f64) -> Result<SisoSosFilter<3, T>, &'static str>
where
    T: Num + Copy + MulAdd<Output = T> + FromPrimitive + ToPrimitive,
{
    // Convert SOS tables from static arrays to slices because new_interpolated expects slices.
    let sos_tables = butter_sos_tables_6::SOS_TABLES
        .each_ref()
        .map(|sec| sec.each_ref().map(|coeffs| &coeffs[..]));
    SisoSosFilter::new_interpolated(
        cutoff_ratio,
        &butter_sos_tables_6::LOG10_CUTOFF_RATIOS,
        sos_tables,
    )
}

// TODO tests

#[cfg(feature = "std")]
#[cfg(test)]
mod test {
    use num_traits::{FromPrimitive, MulAdd, Num, ToPrimitive};
    use crate::sos::SisoSosFilter;
    use super::{butter2, butter4, butter6};
    use super::{butter_sos_tables_2, butter_sos_tables_4, butter_sos_tables_6};

    /// Determine the gain of a filter at a given frequency by simulating its response to a sinewave.
    /// `freq` is normalized to the sample frequency.
    fn simulate_gain_sinewave<const SECTIONS: usize, T>(filter: &mut SisoSosFilter<SECTIONS, T>, freq: f64, n: usize) -> f64 
    where
        T: Num + Copy + MulAdd<Output = T> + FromPrimitive + ToPrimitive,
    {
        let input: Vec<T> = (0..n)
            .map(|i| T::from_f64((2.0 * std::f64::consts::PI * freq * i as f64).sin()).unwrap())
            .collect();
        let original_rms = (input.iter().map(|&v| (v * v).to_f64().unwrap()).sum::<f64>() / n as f64).sqrt();
        let mut output = Vec::with_capacity(n);
        filter.reset();
        for &u in &input {
            output.push(filter.update(u));
        }
        filter.reset();
        let output_rms = (output.iter().map(|&v| (v * v).to_f64().unwrap()).sum::<f64>() / n as f64).sqrt();
        output_rms / original_rms
    }

    fn test_filter<const SECTIONS: usize, T>(
        tname: &str,
        min_cutoff_ratio: f64,
        max_cutoff_ratio: f64,
        butter_fn: fn(f64) -> Result<SisoSosFilter<SECTIONS, T>, &'static str>,
    ) where
        T: Num + Copy + MulAdd<Output = T> + FromPrimitive + ToPrimitive,
    {
        let order = 2 * SECTIONS;
        // Check approximate attenuation at cutoff frequency at the maximum cutoff ratio; should be -3dB or 1/sqrt(2) magnitude
        let mut filtermax = butter_fn(max_cutoff_ratio).unwrap();
        let gain = simulate_gain_sinewave(&mut filtermax, max_cutoff_ratio, 1024);
        let gain_expected = 1.0 / 2f64.sqrt();
        let attenuation_rel_err = (gain - gain_expected).abs() / gain_expected;
        println!("order {order} {tname} attenuation rel err {attenuation_rel_err}");
        assert!(attenuation_rel_err < 0.05);

        // Check convergence of step responses at min and max tabulated cutoff
        let mut filtermin = butter_fn(min_cutoff_ratio).unwrap();
        (0..10_000).for_each(|_| {
            filtermin.update(T::from_f64(1.0).unwrap());
        });
        let step_min_final = filtermin.update(T::from_f64(1.0).unwrap());
        println!("order {order} {tname} step min final {:?}", step_min_final.to_f64().unwrap());
        let step_min_rel_err = (step_min_final.to_f64().unwrap() - 1.0).abs() / 1.0;
        println!("order {order} {tname} step min rel err {step_min_rel_err}");
        assert!(step_min_rel_err < 1e-4);

        filtermax.reset();
        (0..1024).for_each(|_| {
            filtermax.update(T::from_f64(1.0).unwrap());
        });
        let step_max_final = filtermax.update(T::from_f64(1.0).unwrap());
        let step_max_rel_err = (step_max_final.to_f64().unwrap() - 1.0).abs() / 1.0;
        println!("order {order} {tname} step max rel err {step_max_rel_err}");
        assert!(step_max_rel_err < 1e-6);
    }

    #[test]
    fn test_butter2_f32() {
        test_filter::<1, f32>(
            "f32",
            0.005,
            10.0_f64.powf(*butter_sos_tables_2::LOG10_CUTOFF_RATIOS.last().unwrap()),
            butter2,
        );
    }

    #[test]
    fn test_butter2_f64() {
        test_filter::<1, f64>(
            "f64",
            0.002,
            10.0_f64.powf(*butter_sos_tables_2::LOG10_CUTOFF_RATIOS.last().unwrap()),
            butter2,
        );
    }

    #[test]
    fn test_butter4_f32() {
        test_filter::<2, f32>(
            "f32",
            0.005,
            10.0_f64.powf(*butter_sos_tables_4::LOG10_CUTOFF_RATIOS.last().unwrap()),
            butter4,
        );
    }

    #[test]
    fn test_butter4_f64() {
        test_filter::<2, f64>(
            "f64",
            0.0005,
            10.0_f64.powf(*butter_sos_tables_4::LOG10_CUTOFF_RATIOS.last().unwrap()),
            butter4,
        );
    }

    #[test]
    fn test_butter6_f32() {
        test_filter::<3, f32>(
            "f32",
            0.01,
            10.0_f64.powf(*butter_sos_tables_6::LOG10_CUTOFF_RATIOS.last().unwrap()),
            butter6,
        );
    }

    #[test]
    fn test_butter6_f64() {
        test_filter::<3, f64>(
            "f64",
            0.0005,
            10.0_f64.powf(*butter_sos_tables_6::LOG10_CUTOFF_RATIOS.last().unwrap()),
            butter6,
        );
    }
}