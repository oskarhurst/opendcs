/* 
 * Copyright (c) 2018
 * United States Army Corps of Engineers - Hydrologic Engineering Center (USACE/HEC)
 * All Rights Reserved.  USACE PROPRIETARY/CONFIDENTIAL.
 * Source may not be released without written approval from HEC
 */
package decodes.cwms.newresevap;

/**
 *
 * @author RESEVAP program by Steven F. Daly (ERDC/CRREL)
 * conversion to Java by Richard Rachiele (RMA)
 */
public class EvapUtilities
{
    public static final double CONST_K = 0.4;
    public static final double CONST_G = 9.81;

    public EvapUtilities() {
    }

    public static double nu(double airTemp) {
        double nuval = 1.326E-5 * (1.0 + airTemp * (0.006542 + airTemp * (8.301E-6 - 4.84E-9 * airTemp)));
        return nuval;
    }

    public static double spec_ht_air(double airTemp) {
        double spcHeat = 1005.6 + airTemp * (0.017211 + 3.92E-4 * airTemp);
        return spcHeat;
    }

    public static double latent(double airTemp) {
        return airTemp < 0.0 ? (28.34 - 0.00149 * airTemp) * 100000.0 : (25.0 - 0.02274 * airTemp) * 100000.0;
    }

    public static void den_from_rh(double relH, double baroPres, double tempC, double sal, DoubleContainer rhoAdc, DoubleContainer rhoVdc, DoubleContainer q, int iflag) {
        double tK = 273.15;
        double mw = 0.018016;
        double rgas = 8.31441;
        double p0 = 1013.25;
        double t = tempC + 273.15;
        double rhoD = 1.2923 * (273.15 / (tempC + 273.15)) * (baroPres / 1013.25);
        double eSat = satVpr(baroPres, tempC, iflag);
        double e = relH * eSat;
        double fac = 1.0 - 5.37E-4 * sal;
        double rhoV = 100.0 * e * 0.018016 * fac / (8.31441 * t);
        double rhoA = rhoD + rhoV;
        q.d = rhoV / rhoA;
        rhoAdc.d = rhoA;
        rhoVdc.d = rhoV;
    }

    public static double satVpr(double p, double t, int iflag) {
        double[] pa = new double[]{3.46E-6, 4.18E-6};
        double[] pb = new double[]{1.0007, 1.0003};
        double[] e0 = new double[]{6.1121, 6.1115};
        double[] a = new double[]{17.502, 22.452};
        double[] b = new double[]{240.97, 272.55};
        int k = 0;
        if (iflag == 0) {
            if (t < 0.0) {
                k = 1;
            }
        } else if (iflag == 2) {
            k = 1;
        }

        double eSat = (pa[k] * p + pb[k]) * e0[k] * Math.exp(a[k] * t / (b[k] + t));
        return eSat;
    }

    public static double roughness(double r, double cdnr) {
        return r * Math.exp(-0.4 / Math.sqrt(cdnr));
    }

    public static double cdnr(double ru, double z0) {
        double val = Math.log(ru / z0) / 0.4;
        return Math.pow(val, -2.0);
    }

    public static double smooth(double ustar, double gnu) {
        double smoothCoef = 0.135;
        return 0.135 * (gnu / ustar);
    }

    public static double obukhov(double ustar, double tstar, double qstar, double tave, double qave) {
        double tk = 273.15;
        return (tave + 273.15) * ustar * ustar / 3.9240000000000004 / (tstar + 0.61 * (tave + 273.15) * qstar / (1.0 + 0.61 * qave));
    }

    public static double psi_m(double zu, double obukhovLen) {
        double zeta = zu / obukhovLen;
        double term1;
        if (zeta < 0.0) {
            term1 = Math.pow(1.0 - 16.0 * zeta, 0.25);
            return 2.0 * Math.log(0.5 * (1.0 + term1)) + Math.log(0.5 * (1.0 + term1 * term1)) - 2.0 * Math.atan(term1) + 1.570796;
        } else if (zeta == 0.0) {
            return 0.0;
        } else {
            term1 = 0.0;
            if (zeta <= 250.0) {
                term1 = Math.exp(-0.35 * zeta);
            }

            return -(0.7 * zeta + 0.75 * (zeta - 14.3) * term1 + 10.7);
        }
    }

    public static double psi_h(double zs, double obukhovLen) {
        double zeta = zs / obukhovLen;
        double term1;
        if (zeta < 0.0) {
            term1 = Math.pow(1.0 - 16.0 * zeta, 0.25);
            return 2.0 * Math.log(0.5 * (1.0 + term1 * term1));
        } else if (zeta == 0.0) {
            return 0.0;
        } else {
            term1 = 0.0;
            if (zeta <= 250.0) {
                term1 = Math.exp(-0.35 * zeta);
            }

            return -(0.7 * zeta + 0.75 * (zeta - 14.3) * term1 + 10.7);
        }
    }

    public static void coare(double rstar, double gnu, double ustar, DoubleContainer zt, DoubleContainer zq) {
        double[] rs = new double[]{0.135, 0.16, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0};
        double[] at = new double[]{0.177, 1.376, 1.376, 1.026, 1.625, 4.661, 34.904, 1667.19, 588000.0};
        double[] bt = new double[]{0.0, 0.929, 0.929, -0.599, -1.018, -1.475, -2.067, -2.907, -3.935};
        double[] aq = new double[]{0.292, 1.808, 1.808, 1.393, 1.956, 4.994, 30.709, 1448.68, 298000.0};
        double[] bq = new double[]{0.0, 0.826, 0.826, -0.528, -0.87, -1.297, -1.845, -2.682, -3.616};
        double gu = gnu / ustar;
        int idx = 8;

        for(int i = 0; i < 9; ++i) {
            if (rstar < rs[i]) {
                idx = i;
                break;
            }
        }

        zt.d = gu * at[idx] * Math.pow(rstar, bt[idx]);
        zq.d = gu * aq[idx] * Math.pow(rstar, bq[idx]);
    }

    public static double drift_snow(double ustar, double gnu) {
        double z0_smooth = smooth(ustar, gnu);
        double z0_drift = 0.006625891946992864 * ustar * ustar;
        double utmp = (ustar - 0.25) / 0.052;
        double z0_stat = 0.00113 * Math.exp(-(utmp * utmp));
        double z0 = z0_smooth + z0_stat + z0_drift;
        return z0;
    }

    public static void andreas_87(double ustar, double z0, double rstar, DoubleContainer zt, DoubleContainer zq) {
        double[] rs = new double[]{0.135, 2.5, 1000.0};
        double[][] bt = new double[][]{{1.25, 0.0, 0.0}, {0.149, -0.55, 0.0}, {0.317, -0.565, -0.183}};
        double[][] bq = new double[][]{{1.61, 0.0, 0.0}, {0.351, -0.628, 0.0}, {0.396, -0.512, -0.18}};
        int ind_rs = 2;

        for(int i = 0; i < 3; ++i) {
            if (rstar <= rs[i]) {
                ind_rs = i;
                break;
            }
        }

        double lnrstar = Math.log(rstar);
        double lnzsz0 = bt[ind_rs][0] + bt[ind_rs][1] * lnrstar + bt[ind_rs][2] * lnrstar * lnrstar;
        zt.d = z0 * Math.exp(lnzsz0);
        lnzsz0 = bq[ind_rs][0] + bq[ind_rs][1] * lnrstar + bq[ind_rs][2] * lnrstar * lnrstar;
        zq.d = z0 * Math.exp(lnzsz0);
    }

    public static void bulk_coefs(double z0, double zt, double zq, double oblen, double ru, double rt, double rq, DoubleContainer cd, DoubleContainer ch, DoubleContainer ce, int kode) {
        double psim;
        double psiht;
        double psihq;
        if (kode == 0) {
            psim = 0.0;
            psiht = 0.0;
            psihq = 0.0;
        } else {
            psim = psi_m(ru, oblen);
            psiht = psi_h(rt, oblen);
            psihq = psi_h(rq, oblen);
        }

        cd.d = c(ru, ru, z0, z0, psim, psim);
        ch.d = c(ru, rt, z0, zt, psim, psiht);
        ce.d = c(ru, rq, z0, zq, psim, psihq);
    }

    public static double computeFischerUStar(double ur, double rhoAdc, double rhow) {
        double cd;
        if (ur <= 5.0) {
            cd = 0.001;
        } else if (ur >= 15.0) {
            cd = 0.0015;
        } else {
            cd = 0.001 + 5.0E-4 * ((ur - 5.0) / 10.0);
        }

        return Math.sqrt(cd * (rhoAdc / rhow) * ur * ur);
    }

    public static double computeDonelanUStar(double ur, double rhoAdc, double rhow) {
        double cd = (0.37 + 0.137 * ur) / 1000.0;
        return Math.sqrt(cd * rhoAdc / rhow) * ur;
    }

    public static double c(double z1, double z2, double z01, double z02, double psi1, double psi2) {
        double cval = 0.16000000000000003 / ((Math.log(z1 / z01) - psi1) * (Math.log(z2 / z02) - psi2));
        return cval;
    }

    public static double speed(double ur, double ustar, double zi, double l) {
        double w0 = 0.5;
        double beta = 1.25;
        if (Math.abs(l) >= 1000.0) {
            return ur;
        } else {
            double spd;
            if (l < 0.0) {
                spd = ustar * Math.pow(-zi / (0.4 * l), 0.333333);
                double spd = Math.sqrt(ur * ur + 1.25 * spd * 1.25 * spd);
                return spd;
            } else {
                spd = ur + 0.5;
                return spd;
            }
        }
    }

    public static void fluxes(double ustar, double tstar, double qstar, double rhoa, double cp, double lv, DoubleContainer tau, DoubleContainer hs, DoubleContainer hl, DoubleContainer evap) {
        double rhow = 1000.0;
        tau.d = rhoa * ustar * ustar;
        hs.d = -rhoa * cp * ustar * tstar;
        hl.d = -rhoa * lv * ustar * qstar;
        evap.d = hl.d / (lv * 1000.0);
        evap.d = 8.64E7 * evap.d;
    }

    public static double den_h2o(double tx) {
        double den = 1000.0 - 0.019549 * Math.pow(Math.abs(tx - 4.0), 1.68);
        return den;
    }

    public static double cp_h2o(double px) {
        double tx;
        if (px < 0.0) {
            tx = 0.0;
        } else {
            tx = px;
        }

        double rx = 34.5 - tx;
        double cp = 4174.9 + 1.6659 * (Math.exp(rx / 10.6) + Math.exp(-1.0 * rx / 10.6));
        return cp;
    }

    public static double f_bo_star(double tempC, double baroPres, double sal, int iflag) {
        double tK = 273.15;
        double[] a = new double[]{17.502, 22.452};
        double[] b = new double[]{240.97, 272.55};
        double lv = latent(tempC);
        double cp = spec_ht_air(tempC);
        double f = 1.0;
        DoubleContainer rhoA = new DoubleContainer();
        DoubleContainer rhoVS = new DoubleContainer();
        DoubleContainer dum = new DoubleContainer();
        den_from_rh(f, baroPres, tempC, sal, rhoA, rhoVS, dum, iflag);
        double rhoD = rhoA.d + rhoVS.d;
        double bo_fac = a[iflag - 1] * b[iflag - 1] / Math.pow(b[iflag - 1] + tempC, 2.0) - 1.0 / (tempC + 273.15);
        return rhoD * cp / (lv * rhoVS.d * bo_fac);
    }

    public static double cdd(double x) {
        return (0.37 + 0.137 * x) * 0.001;
    }

    public static class DoubleContainer {
        double d;

        public DoubleContainer() {
            this.d = 0.0;
        }

        public DoubleContainer(double d) {
            this.d = d;
        }
    }
    
}
