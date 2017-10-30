import sys,os
import numpy as np
from scipy import interpolate
from scipy import stats

def load_frame(file_name):
    measurements = []
    header = []
    frame_file = open(file_name, "r")
    content = frame_file.readlines()
    for line in content:
        if line.startswith("#"):
            header.append(line)
        numbers_str = line.split()
        try:
             numbers = [float(x) for x in numbers_str]
        except:
            header.append(line)
            continue
        if len(numbers) > 0:
            measurements.append(numbers)
    frame_file.close()
    return measurements, header

def calc_chi2(curve_exp,sigma_exp, curve_calc):
    chi2 = np.sum(np.square(np.divide((curve_calc-curve_exp),sigma_exp)))/(len(curve_calc)-1)
    return chi2

def fit2(curve_exp,sigma_exp, curve_calc):
    s2 = np.square(sigma_exp)
    sum1 = np.divide(np.multiply(curve_exp,curve_calc),s2)
    sum2 = np.divide(np.square(curve_calc),s2)
    c = np.sum(sum1)/np.sum(sum2)
    return curve_calc*c

def calc_chi_distr2(curve_exp,sigma_exp,int_curves):
    chi2_vals = []
    for calc_curve in int_curves:
        ls_fit2 = fit2(curve_exp,sigma_exp,calc_curve)
        chi2_val = calc_chi2(curve_exp,sigma_exp,ls_fit2)
        chi2_vals.append(chi2_val)
    return np.array(chi2_vals)


def run(args):
    if len(args) < 1:
        print "\nUsage: python calc_z.py <data.dat> [model.pdb]"
        sys.exit()

    data_dir = os.path.dirname(os.path.realpath(__file__))
    loaded = np.load(os.path.join(data_dir, 'curves.npz'))
    I_calc = loaded['curves']
    s_calc = loaded['range']

    filename_data = args[0]

    filename_model = ''
    if len(args) > 1:
        filename_model = args[1]

    frame, header= load_frame(filename_data)
    np_data = np.array(frame)
    s_exp = np_data[:,0]
    I_exp = np_data[:,1]
    I_sigma = np_data[:,2]

    if max(s_exp) > 2:
        s_exp /=10

    if max(s_exp) > 0.5:
        ind_remain = s_exp<0.5
        s_exp = s_exp[ind_remain]
        I_exp = I_exp[ind_remain]
        I_sigma = I_sigma[ind_remain]

    eps = 0.000001
    f = interpolate.interp1d(s_calc, I_calc, kind='linear')
    int_curves = f(s_exp)
    chi_distribution = calc_chi_distr2(I_exp, I_sigma, int_curves)

    chi_max = np.max(chi_distribution)+eps
    zero_level = 1 / chi_max
    chi_distribution = chi_distribution / chi_max
    parameters = stats.beta.fit(chi_distribution, floc=0, fscale=1)
    distr_frozen = stats.beta.freeze(*parameters)
    distr_cdf = distr_frozen.cdf(zero_level)
    z_limit = -stats.norm(0, 1).ppf(distr_cdf)
    print "Data Z0 = %0.2f" % z_limit

    if filename_model:
        try:
            import IMP
            import IMP.atom
            import IMP.core
            import IMP.saxs
            import IMP.foxs
        except:
            print "Model-data fitting requires IMP."
            sys.exit()

        m = IMP.Model()
        mps = IMP.atom.read_multimodel_pdb(filename_model, m, IMP.atom.NonWaterNonHydrogenPDBSelector(), True)
        particles = []
        for mp in mps:
            particles += IMP.atom.get_by_type(mp, IMP.atom.ATOM_TYPE)

        ft = IMP.saxs.get_default_form_factor_table()
        for i in range(0, len(particles)):
            radius = ft.get_radius(particles[i])
            IMP.core.XYZR.setup_particle(particles[i], radius)

        s = IMP.saxs.SolventAccessibleSurface()
        surface_area = s.get_solvent_accessibility(IMP.core.XYZRs(particles))

        delta_q = 0.5 / 1000
        model_profile = IMP.saxs.Profile(0.0, 0.5, delta_q)
        model_profile.calculate_profile_partial(particles, surface_area)

        model_qs = model_profile.get_qs()
        model_ints = model_profile.get_intensities()
        model_ints = np.array(model_ints)/np.median(model_ints)
        f = interpolate.interp1d(model_qs, model_ints, kind='quadratic')
        curve_calc = f(s_exp)

        ls_fit = fit2(I_exp, I_sigma, curve_calc)
        chi2 = calc_chi2(I_exp, I_sigma, ls_fit)
        chi_tranformed = chi2 / chi_max
        distr_cdf = distr_frozen.cdf(chi_tranformed)
        z = -stats.norm(0, 1).ppf(distr_cdf)
        print "Model chi2 = %0.2f" % chi2
        print "Model-data Z-score = %0.2f" % z

if (__name__ == "__main__"):
    run(args=sys.argv[1:])





