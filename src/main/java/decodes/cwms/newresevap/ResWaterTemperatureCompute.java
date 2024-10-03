/* 
 * Copyright (c) 2018
 * United States Army Corps of Engineers - Hydrologic Engineering Center (USACE/HEC)
 * All Rights Reserved.  USACE PROPRIETARY/CONFIDENTIAL.
 * Source may not be released without written approval from HEC
 */
package decodes.cwms.newresevap;

//import hec.heclib.util.HecTime;

import decodes.cwms.resevapcalc.Const;
import decodes.cwms.resevapcalc.EvapReservoir;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Date;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *  Class used to compute reservoir temperature profile
 */
public class ResWaterTemperatureCompute
{
    private static final Logger LOGGER = Logger.getLogger(ResWtCompute.class.getName());
    BufferedWriter _tout = null;
    EvapReservoir _reservoir;
    int _resj_old = -1;
    private double[] a;
    private double[] b;
    private double[] c;
    private double[] r;
    private double[] u;
    private double[] rhow_org;
    private double[] cp_org;
    private double[] wt_tmp;
    private double[] pe_mix_out;
    final int NMAX = 500;
    double[] gam = new double[500];

    public ResWtCompute(EvapReservoir reservoir) {
        this._reservoir = reservoir;
        this.a = new double[1000];
        this.b = new double[1000];
        this.c = new double[1000];
        this.r = new double[1000];
        this.u = new double[1000];
        this.rhow_org = new double[1000];
        this.cp_org = new double[1000];
        this.wt_tmp = new double[1000];
        this.pe_mix_out = new double[1000];
    }

    public void setOutfile(BufferedWriter tout) {
        this._tout = tout;
    }

    public boolean computeReservoirTemp(HecTime currentTime, MetComputation metComputation, double delT) {
        Arrays.fill(this.a, 0.0);
        Arrays.fill(this.b, 0.0);
        Arrays.fill(this.c, 0.0);
        Arrays.fill(this.r, 0.0);
        Arrays.fill(this.u, 0.0);
        Arrays.fill(this.rhow_org, 0.0);
        Arrays.fill(this.cp_org, 0.0);
        Arrays.fill(this.wt_tmp, 0.0);
        Arrays.fill(this.pe_mix_out, 0.0);
        double[] zd = this._reservoir._zd;
        double[] kz = this._reservoir._kz;
        double[] zarea = this._reservoir._zarea;
        double[] delz = this._reservoir._delz;
        double[] zvol = this._reservoir._zvol;
        double[] ztop = this._reservoir._ztop;
        double[] rhow = this._reservoir._rhow;
        double[] cp = this._reservoir._cp;
        double[] wt = this._reservoir._wt;
        WindShearMethod windShearMethod = this._reservoir._windShearMethod;
        double thermalDiffusivityCoefficient = this._reservoir._thermalDiffusivityCoefficient;
        double wsel = this._reservoir.getElevation();
        double surfArea = this._reservoir._surfArea;
        double grav = 9.81;
        double SOLAR = metComputation._solar;
        double flxir = metComputation._flxir;
        double flxir_out = metComputation._flxir_out;
        double hs = metComputation._evapWater._hs;
        double hl = metComputation._evapWater._hl;
        double evap = metComputation._evapWater._evap;
        double ustar = metComputation._evapWater._ustar;
        double katten = this._reservoir._attenuationConst;
        double ur = metComputation._metData._windSpeed_current;
        double rh = metComputation._metData._relHumidity_current;
        double tr = metComputation._metData._airTemp_current;
        double p = metComputation._metData._airPressure_current;
        double theta = 1.0;
        double albedo = 0.08;
        double penfrac = 0.4;
        double eta_conv = 0.5;
        double wind_critic = 1.0;
        double eta_stir = 0.4;
        int resj_old = this._reservoir._resj_old;
        int resj = this._reservoir.getResj();
        double avg_rhow = 0.0;
        double total_vol = 0.0;
        double sum_energy_in = 0.0;
        double TKE = 0.0;
        int i;
        if (resj != resj_old) {
            for(i = 0; i < resj_old; ++i) {
                this.wt_tmp[i] = wt[i];
            }

            i = resj;

            for(int j = resj_old; i >= 0 && j >= 0; --j) {
                wt[i] = this.wt_tmp[j];
                --i;
            }

            if (resj > resj_old) {
                for(int k = i; k >= 0; --k) {
                    wt[k] = this.wt_tmp[0];
                }
            }
        }

        for(i = 0; i <= resj; ++i) {
            if (wt[i] < 0.0) {
                wt[i] = 0.0;
            }

            rhow[i] = EvapUtilities.den_h2o(wt[i]);
            cp[i] = EvapUtilities.cp_h2o(wt[i]);
            avg_rhow += rhow[i] * zvol[i];
            total_vol += zvol[i];
            sum_energy_in += rhow[i] * cp[i] * zvol[i] * wt[i];
        }

        avg_rhow /= total_vol;

        for(i = resj; i >= 1; --i) {
            double sFreq;
            if (rhow[i] - rhow[i - 1] != 0.0) {
                sFreq = grav / avg_rhow * Math.abs(rhow[i] - rhow[i - 1]) / (zd[i - 1] - zd[i]);
            } else {
                sFreq = 7.0E-5;
            }

            if (sFreq < 7.0E-5) {
                sFreq = 7.0E-5;
            }

            double surfArea_x;
            if (surfArea < 350.0) {
                surfArea_x = surfArea;
            } else {
                surfArea_x = 350.0;
            }

            double kzx = 8.17E-4 * Math.pow(surfArea_x, 0.56) * Math.pow(Math.abs(sFreq), -0.43);
            kzx = 1.0E-4 * kzx;
            kz[i] = thermalDiffusivityCoefficient * kzx * cp[i] * rhow[i];
        }

        for(i = 0; i <= resj; ++i) {
            double al;
            double abar;
            double au;
            if (i == 0) {
                au = zarea[i];
                abar = 0.5 * zarea[i];
                al = 0.0;
            } else {
                au = zarea[i];
                abar = 0.5 * (zarea[i] + zarea[i - 1]);
                al = zarea[i - 1];
            }

            double delZu;
            double dhatZu;
            if (i < resj) {
                delZu = 0.5 * (delz[i + 1] + delz[i]);
                dhatZu = (kz[i + 1] / (cp[i + 1] * rhow[i + 1]) * delz[i + 1] + kz[i] / (cp[i] * rhow[i]) * delz[i]) / (delz[i + 1] + delz[i]);
            } else {
                dhatZu = 0.0;
                delZu = 1.0;
            }

            double delZl;
            double dhatZl;
            if (i > 0) {
                delZl = 0.5 * (delz[i - 1] + delz[i]);
                dhatZl = (kz[i - 1] / (cp[i - 1] * rhow[i - 1]) * delz[i - 1] + kz[i] / (cp[i] * rhow[i]) * delz[i]) / (delz[i - 1] + delz[i]);
            } else {
                dhatZl = 0.0;
                delZl = 1.0;
            }

            if (i > 0) {
                this.a[i] = -1.0 * delT / delz[i] * dhatZl / delZl * theta * al / abar;
            } else {
                this.a[i] = 0.0;
            }

            this.b[i] = 1.0 + delT / delz[i] * dhatZu / delZu * theta * au / abar + delT / delz[i] * dhatZl / delZl * theta * al / abar;
            if (i < resj) {
                this.c[i] = -1.0 * delT / delz[i] * dhatZu / delZu * theta * au / abar;
            } else {
                this.c[i] = 0.0;
            }

            if (i == 0) {
                this.r[i] = wt[i] + delT / delz[i] * au / abar * dhatZu / delZu * (1.0 - theta) * (wt[i + 1] - wt[i]);
            }

            if (i > 0 && i < resj) {
                this.r[i] = wt[i] + delT / delz[i] * au / abar * (dhatZu / delZu * (1.0 - theta) * (wt[i + 1] - wt[i]) - dhatZl / delZl * al / abar * (1.0 - theta) * (wt[i] - wt[i - 1]));
            }

            if (i == resj) {
                this.r[i] = wt[i] + delT / delz[i] * al / abar * dhatZl / delZl * (1.0 - theta) * (wt[i] - wt[i - 1]);
            }
        }

        double topDepth = 0.0;
        double solar_tot = 0.0;

        double fi;
        for(i = resj; i > 1; --i) {
            double bottomDepth = topDepth + delz[i];
            fi = SOLAR * penfrac * (1.0 - albedo) * (Math.exp(-1.0 * katten * topDepth) * zarea[i] - Math.exp(-1.0 * katten * bottomDepth) * zarea[i - 1]) / zvol[i];
            solar_tot += SOLAR * penfrac * (1.0 - albedo) * (Math.exp(-1.0 * katten * topDepth) * zarea[i] - Math.exp(-1.0 * katten * bottomDepth) * zarea[i - 1]) / zvol[i];
            this.r[i] += fi * delT / (cp[i] * rhow[i]);
            topDepth = bottomDepth;
        }

        fi = (SOLAR * (1.0 - penfrac) * (1.0 - albedo) + flxir + flxir_out - hs - hl) * zarea[resj] / zvol[resj];
        double var10000 = solar_tot + SOLAR * (1.0 - penfrac) * (1.0 - albedo) * zarea[resj];
        this.r[resj] += fi * delT / (cp[resj] * rhow[resj]);
        double fi_check = (SOLAR * (1.0 - albedo) + flxir + flxir_out - hs - hl) * delT;
        fi_check *= zarea[resj];
        this.tridag(this.a, this.b, this.c, this.r, this.u, resj + 1);
        double tempDepth = 0.25;
        double sum_energy_diff = 0.0;

        for(i = resj; i >= 0; --i) {
            wt[i] = this.u[i];
            rhow[i] = EvapUtilities.den_h2o(wt[i]);
            this.rhow_org[i] = rhow[i];
            cp[i] = EvapUtilities.cp_h2o(wt[i]);
            this.cp_org[i] = cp[i];
            sum_energy_diff += wt[i] * cp[i] * rhow[i] * zvol[i];
            if (i > 0) {
                tempDepth += 0.5 * (delz[i] + delz[i - 1]);
            }
        }

        EvapUtilities.DoubleContainer rhoAdc = new EvapUtilities.DoubleContainer();
        EvapUtilities.DoubleContainer rhoVdc = new EvapUtilities.DoubleContainer();
        EvapUtilities.DoubleContainer spHumdc = new EvapUtilities.DoubleContainer();
        int sml_flag = true;
        int sml = resj;
        double sml_delz = delz[resj];
        double sml_vol = zvol[resj];
        int i = resj;

        for(int iend = true; i >= 2 && iend; --i) {
            rhow[i] = EvapUtilities.den_h2o(wt[i]);
            cp[i] = EvapUtilities.cp_h2o(wt[i]);
            double wt_mix = (wt[i] * sml_vol * cp[i] + wt[i - 1] * zvol[i - 1] * cp[i - 1]) / (sml_vol * cp[i] + zvol[i - 1] * cp[i - 1]);
            double rhow_mix = EvapUtilities.den_h2o(wt_mix);
            double z_mix = this.zcom(resj, sml - 1, rhow_mix);
            double z_sml = this.zcom(resj, sml, -1.0);
            double z_next = this.zcom(sml - 1, sml - 1, -1.0);
            double pe_mix = grav * (rhow_mix * (sml_vol + zvol[i - 1]) * (z_mix - ztop[i - 2]) - (rhow[i] * sml_vol * (z_sml - ztop[i - 2]) + rhow[i - 1] * zvol[i - 1] * (z_next - ztop[i - 2])));
            this.pe_mix_out[i] = pe_mix;
            int j;
            if (pe_mix < 0.0) {
                sml = i - 1;
                sml_delz += delz[i - 1];
                sml_vol += zvol[i - 1];

                for(j = i - 1; j <= resj; ++j) {
                    wt[j] = wt_mix;
                }
            } else {
                if (sml_flag) {
                    sml_flag = false;
                    double org_PE = 0.0;
                    double new_PE = 0.0;

                    for(j = sml; j <= resj; ++j) {
                        org_PE += this.rhow_org[j] * delz[j] * (ztop[j] + ztop[i - 1]) / 2.0;
                        new_PE += this.rhow_org[j] * delz[j];
                    }

                    new_PE = new_PE * (ztop[resj] + ztop[sml - 1]) / 2.0;
                    double w_3 = grav / (rhow[resj] * delT) * (org_PE - new_PE);
                    if (w_3 < 0.0) {
                        w_3 = 0.0;
                    }

                    double KE_conv = eta_conv * rhow[resj] * zarea[sml - 1] * w_3 * delT;
                    double u_H2O_star = 0.0;
                    double KE_stir;
                    if (ur > wind_critic) {
                        int IFLAG = 1;
                        double Salinity = 0.0;
                        EvapUtilities.den_from_rh(rh, p, tr, Salinity, rhoAdc, rhoVdc, spHumdc, IFLAG);
                        if (WindShearMethod.DONELAN.equals(windShearMethod)) {
                            u_H2O_star = EvapUtilities.computeDonelanUStar(ur, rhoAdc.d, rhow[resj]);
                        } else {
                            u_H2O_star = EvapUtilities.computeFischerUStar(ur, rhoAdc.d, rhow[resj]);
                        }

                        KE_stir = eta_stir * rhow[resj] * zarea[resj] * Math.pow(u_H2O_star, 3.0) * delT;
                    } else {
                        KE_stir = 0.0;
                    }

                    if (KE_stir < 0.0) {
                        String msg = "KE_stir < 0.0 \neta_stir =" + eta_stir + "\nrhow(resj) =" + rhow[resj] + "\nzarea(resj) =" + zarea[resj] + "\nu_H2O_star =" + u_H2O_star + "\ndelT =" + delT;
                        Logger.getLogger(ResWtCompute.class.getName()).log(Level.SEVERE, msg);
                        return false;
                    }

                    TKE = KE_stir + KE_conv;
                }

                if (TKE >= pe_mix) {
                    sml = i - 1;
                    sml_delz += delz[i - 1];
                    sml_vol += zvol[i - 1];

                    for(j = i - 1; j <= resj; ++j) {
                        wt[j] = wt_mix;
                    }

                    TKE -= pe_mix;
                } else {
                    iend = false;
                }
            }
        }

        double sum_energy_mix = 0.0;

        for(i = resj; i >= 0; --i) {
            rhow[i] = EvapUtilities.den_h2o(wt[i]);
            cp[i] = EvapUtilities.cp_h2o(wt[i]);
            sum_energy_mix += wt[i] * cp[i] * rhow[i] * zvol[i];
        }

        var10000 = wt[resj];
        double changeEng = sum_energy_diff - sum_energy_in;
        double engBalance = changeEng - fi_check;
        double efficiency = engBalance / fi_check * 100.0;
        String dateTimeStr = currentTime.date(4) + "     " + currentTime.getTime(false);
        NumberFormat nfi_4 = new NumberFormat("%4d");
        NumberFormat nf9_2 = new NumberFormat("%9.2f");
        NumberFormat nfe14_8 = new NumberFormat("%14.8E");
        String vals1 = nfi_4.form((long)(sml + 1)) + " " + nfi_4.form((long)(resj + 1)) + " " + nf9_2.form(wsel) + " ";
        String vals2 = nfe14_8.form(sum_energy_in) + " " + nfe14_8.form(sum_energy_diff) + " " + nfe14_8.form(fi_check) + " " + nfe14_8.form(changeEng) + " " + nfe14_8.form(engBalance) + " " + nfe14_8.form(efficiency) + " " + nfe14_8.form(zarea[resj]) + " ";
        if (this._tout != null) {
            try {
                this._tout.write(dateTimeStr + vals1 + vals2);
                this._tout.newLine();
            } catch (IOException var168) {
                IOException ioe = var168;
                StringBuilder msg = new StringBuilder(dateTimeStr + System.lineSeparator());
                msg.append("resj ").append(resj).append(System.lineSeparator());
                msg.append("sum_energy_in ").append(sum_energy_in).append(System.lineSeparator());
                msg.append("sum_energy_diff ").append(sum_energy_diff).append(System.lineSeparator());
                msg.append("fi_check ").append(fi_check).append(System.lineSeparator());
                msg.append("changeEng ").append(changeEng).append(System.lineSeparator());
                msg.append("engBalance ").append(engBalance).append(System.lineSeparator());
                msg.append("efficiency ").append(efficiency).append(System.lineSeparator());
                msg.append("zarea[resj] ").append(zarea[resj]).append(System.lineSeparator());

                for(i = resj; i >= 0; --i) {
                    msg.append(" i, wt[i], zvol[i] ").append(wt[i]).append("  ").append(zvol[i]).append(System.lineSeparator());
                }

                LOGGER.log(Level.SEVERE, ioe, msg::toString);
                return false;
            }
        }

        return true;
    }

    protected double zcom(int itop, int ibottom, double xrhow) {
        double[] rhow = this._reservoir._rhow;
        double total = 0.0;
        double totalpvol = 0.0;
        double[] zvol = this._reservoir._zvol;
        double[] ztop = this._reservoir._ztop;

        for(int j = ibottom; j <= itop; ++j) {
            double zrhow;
            if (xrhow == -1.0) {
                zrhow = rhow[j];
            } else {
                zrhow = xrhow;
            }

            if (j > 0) {
                total += zrhow * zvol[j] * (ztop[j] + ztop[j - 1]) / 2.0;
            } else {
                total += zrhow * zvol[j] * (ztop[j] + 0.0) / 2.0;
            }

            totalpvol += zrhow * zvol[j];
        }

        if (totalpvol != 0.0) {
            double zout = total / totalpvol;
            return zout;
        } else {
            return 0.0;
        }
    }

    private boolean tridag(double[] a, double[] b, double[] c, double[] r, double[] u, int n) {
        if (b[0] == 0.0) {
            Logger.getLogger(ResWtCompute.class.getName()).log(Level.SEVERE, "tridag: rewrite equations");
            return false;
        } else {
            double bet = b[0];
            u[0] = r[0] / bet;

            int j;
            for(j = 1; j < n; ++j) {
                this.gam[j] = c[j - 1] / bet;
                bet = b[j] - a[j] * this.gam[j];
                if (bet == 0.0) {
                    Logger.getLogger(ResWtCompute.class.getName()).log(Level.SEVERE, " tridag failed");
                    return false;
                }

                u[j] = (r[j] - a[j] * u[j - 1]) / bet;
            }

            for(j = n - 2; j >= 0; --j) {
                u[j] -= this.gam[j + 1] * u[j + 1];
            }

            return true;
        }
    }
    
    
}
