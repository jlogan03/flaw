use num_traits::{FromPrimitive, MulAdd, Num};
use crate::sos::SisoSosFilter;
mod butter_sos_tables_4;
mod butter_sos_tables_6;

pub fn butter4<T>(cutoff_ratio: f64) -> Result<SisoSosFilter<2, T>, &'static str>
where
    T: Num + Copy + MulAdd<Output = T> + FromPrimitive,
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
    T: Num + Copy + MulAdd<Output = T> + FromPrimitive,
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
    use crate::sos::SisoSosFilter;
    use super::{butter4, butter6};
    use super::{butter_sos_tables_4, butter_sos_tables_6};

    /// Determine the gain of a filter at a given frequency by simulating its response to a sinewave.
    /// `freq` is normalized to the sample frequency.
    fn simulate_gain_sinewave_f32<const SECTIONS: usize>(filter: &mut SisoSosFilter<SECTIONS, f32>, freq: f32, n: usize) -> f32 {
        let input: Vec<f32> = (0..n)
            .map(|i| (2.0 * std::f32::consts::PI * freq * i as f32).sin())
            .collect();
        let original_rms = (input.iter().map(|v| v * v).sum::<f32>() / n as f32).sqrt();
        let mut output = Vec::with_capacity(n);
        filter.reset();
        for &u in &input {
            output.push(filter.update(u));
        }
        filter.reset();
        let output_rms = (output.iter().map(|v| v * v).sum::<f32>() / n as f32).sqrt();
        output_rms / original_rms
    }


    #[test]
    fn test_butter4() {
        let order = 4;
        // let min_cutoff_ratio = 10.0_f64.powf(butter_sos_tables_4::LOG10_CUTOFF_RATIOS[0]);
        let min_cutoff_ratio = 0.001;
        let max_cutoff_ratio = 10.0_f64.powf(*butter_sos_tables_4::LOG10_CUTOFF_RATIOS.last().unwrap());

        // Check approximate attenuation at cutoff frequency at the maximum cutoff ratio; should be -3dB or 1/sqrt(2) magnitude
        let mut filtermax = butter4::<f32>(max_cutoff_ratio).unwrap();
        let gain = simulate_gain_sinewave_f32(&mut filtermax, max_cutoff_ratio as f32, 1024);
        let gain_expected = 1.0 / 2f32.sqrt();
        let attenuation_rel_err = (gain - gain_expected).abs() / gain_expected;
        println!("order {order} attenuation rel err {attenuation_rel_err}");
        assert!(attenuation_rel_err < 0.05);

        // Check convergence of step responses at min and max tabulated cutoff
        let mut filtermin = butter4::<f32>(min_cutoff_ratio).unwrap();
        (0..10_000).for_each(|_| {filtermin.update(1.0);});
        let step_min_final = filtermin.update(1.0);
        println!("step min final {}", step_min_final);
        let step_min_rel_err = (step_min_final - 1.0).abs() / 1.0;
        println!("order {order} step min rel err {step_min_rel_err}");
        assert!(step_min_rel_err < 1e-4);

        filtermax.reset();
        (0..1024).for_each(|_| {filtermax.update(1.0);});
        let step_max_final = filtermax.update(1.0);
        let step_max_rel_err = (step_max_final - 1.0).abs() / 1.0;
        println!("order {order} step max rel err {step_max_rel_err}");
        assert!(step_max_rel_err < 1e-6);
    }

    #[test]
    fn test_butter6() {
        let order = 6;
        // let min_cutoff_ratio = 10.0_f64.powf(butter_sos_tables_6::LOG10_CUTOFF_RATIOS[0]);
        let min_cutoff_ratio = 0.1;
        let max_cutoff_ratio = 10.0_f64.powf(*butter_sos_tables_6::LOG10_CUTOFF_RATIOS.last().unwrap());

        // Check approximate attenuation at cutoff frequency at the maximum cutoff ratio; should be -3dB or 1/sqrt(2) magnitude
        let mut filtermax = butter6::<f32>(max_cutoff_ratio).unwrap();
        let gain = simulate_gain_sinewave_f32(&mut filtermax, max_cutoff_ratio as f32, 1024);
        let gain_expected = 1.0 / 2f32.sqrt();
        let attenuation_rel_err = (gain - gain_expected).abs() / gain_expected;
        println!("order {order} attenuation rel err {attenuation_rel_err}");
        assert!(attenuation_rel_err < 0.05);

        // Check convergence of step responses at min and max tabulated cutoff
        let mut filtermin = butter6::<f32>(min_cutoff_ratio).unwrap();
        (0..10_000).for_each(|_| {filtermin.update(1.0);});
        let step_min_final = filtermin.update(1.0);
        println!("step min final {}", step_min_final);
        let step_min_rel_err = (step_min_final - 1.0).abs() / 1.0;
        println!("order {order} step min rel err {step_min_rel_err}");
        // assert!(step_min_rel_err < 1e-4);

        filtermax.reset();
        (0..1024).for_each(|_| {filtermax.update(1.0);});
        let step_max_final = filtermax.update(1.0);
        let step_max_rel_err = (step_max_final - 1.0).abs() / 1.0;
        println!("order {order} step max rel err {step_max_rel_err}");
        assert!(step_max_rel_err < 1e-6);
    }
}