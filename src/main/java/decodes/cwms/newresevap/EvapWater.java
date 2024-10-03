/* 
 * Copyright (c) 2018
 * United States Army Corps of Engineers - Hydrologic Engineering Center (USACE/HEC)
 * All Rights Reserved.  USACE PROPRIETARY/CONFIDENTIAL.
 * Source may not be released without written approval from HEC
 */
package decodes.cwms.newresevap;

import decodes.cwms.resevapcalc.Const;

/**
 *
 * @author RESEVAP program by Steven F. Daly  (ERDC/CRREL)
 * conversion to Java by Richard Rachiele (RMA)
 * 
 * Class used to compute and hold evaporation rate, latent and sensible heat
 * fluxes.
 */
public class EvapWater
{
    double _evap;
    double _hs;
    double _hl;
    double _ustar;
    double _tstar;
    double _pstar;
    double _qstar;
    double _rstar;
    double _obukhovLen;

    public EvapWater() {
    }

    public void evap_water(double surfaceTemp, double airTemp, double relHumidity, double windspeed, double airPressure, double rt, double ru, double rq) {
        double s = 0.0;
        int iflag = 1;
        double r10 = 10.0;
        double zi = 600.0;
        double conv = 0.001;
        int ice = true;
        double err_lim = -900.0;
        double ur = windspeed;
        double tr = airTemp;
        double rh = relHumidity;
        double p = airPressure;
        double ts = surfaceTemp;
        if (!(ur < -900.0) && !(tr < -900.0) && !(rh < -900.0) && !(p < -900.0)) {
            double thetar = tr + 9.81 / EvapUtilities.spec_ht_air(tr) * rt;
            if (ts < 0.0) {
                ts = 0.0;
            }

            if (ur < 0.1) {
                ur = 0.1;
            }

            double del_theta = ts - thetar;
            double tave = 0.5 * (ts + tr);
            double rhs = 1.0;
            EvapUtilities.DoubleContainer rhoAdc = new EvapUtilities.DoubleContainer();
            EvapUtilities.DoubleContainer qsdc = new EvapUtilities.DoubleContainer();
            EvapUtilities.DoubleContainer qrdc = new EvapUtilities.DoubleContainer();
            EvapUtilities.DoubleContainer dumdc = new EvapUtilities.DoubleContainer();
            EvapUtilities.DoubleContainer dum1dc = new EvapUtilities.DoubleContainer();
            EvapUtilities.den_from_rh(rhs, p, ts, 0.0, rhoAdc, dumdc, qsdc, iflag);
            double qs = qsdc.d;
            double rhoa = rhoAdc.d;
            double f = rh / 100.0;
            EvapUtilities.den_from_rh(f, p, tr, 0.0, dumdc, dum1dc, qrdc, iflag);
            double qr = qrdc.d;
            double del_q = qs - qr;
            double qave = 0.5 * (qs + qr);
            double gnu = EvapUtilities.nu(ts);
            int its = true;
            double cdnr = EvapUtilities.cdd(ur);
            double ustar = ur * Math.sqrt(cdnr);
            if (ustar < 0.01) {
                ustar = 0.01;
            }

            double z0 = EvapUtilities.roughness(ru, cdnr) + EvapUtilities.smooth(ustar, gnu);
            EvapUtilities.DoubleContainer ztdc = new EvapUtilities.DoubleContainer();
            EvapUtilities.DoubleContainer zqdc = new EvapUtilities.DoubleContainer();
            double rstar = ustar * z0 / gnu;
            EvapUtilities.coare(rstar, gnu, ustar, ztdc, zqdc);
            double zt = ztdc.d;
            double zq = zqdc.d;
            EvapUtilities.DoubleContainer cddc = new EvapUtilities.DoubleContainer();
            EvapUtilities.DoubleContainer cedc = new EvapUtilities.DoubleContainer();
            EvapUtilities.DoubleContainer chdc = new EvapUtilities.DoubleContainer();
            int kode = 0;
            double dum = 0.0;
            EvapUtilities.bulk_coefs(z0, zt, zq, dum, ru, rt, rq, cddc, chdc, cedc, kode);
            double cd = cddc.d;
            double ch = chdc.d;
            double ce = cedc.d;
            ustar = ur * Math.sqrt(cd);
            if (ustar < 0.01) {
                ustar = 0.01;
            }

            double tstar = -(ch * ur * del_theta) / ustar;
            double qstar = -(ce * ur * del_q) / ustar;
            double l = EvapUtilities.obukhov(ustar, tstar, qstar, tave, qave);
            kode = 2;
            int maxiter = 20;

            for(int its = 2; its <= maxiter; ++its) {
                double ustar_old = ustar;
                double tstar_old = tstar;
                double qstar_old = qstar;
                EvapUtilities.bulk_coefs(z0, zt, zq, l, 10.0, rt, rq, cddc, chdc, cedc, kode);
                cd = cddc.d;
                ch = chdc.d;
                ce = cedc.d;
                double u10 = ustar / Math.sqrt(cd);
                double cdn10 = EvapUtilities.cdd(u10);
                z0 = EvapUtilities.roughness(10.0, cdn10) + EvapUtilities.smooth(ustar, gnu);
                rstar = ustar * z0 / gnu;
                EvapUtilities.coare(rstar, gnu, ustar, ztdc, zqdc);
                zt = ztdc.d;
                zq = zqdc.d;
                EvapUtilities.bulk_coefs(z0, zt, zq, l, ru, rt, rq, cddc, chdc, cedc, kode);
                cd = cddc.d;
                ch = chdc.d;
                ce = cedc.d;
                double wind = EvapUtilities.speed(ur, ustar, 600.0, l);
                ustar = wind * Math.sqrt(cd);
                if (ustar < 0.01) {
                    ustar = 0.01;
                }

                tstar = -(ch * wind * del_theta) / ustar;
                qstar = -(ce * wind * del_q) / ustar;
                l = EvapUtilities.obukhov(ustar, tstar, qstar, tave, qave);
                double testu = Math.abs((ustar - ustar_old) / ustar);
                double testt;
                if (tstar != 0.0) {
                    testt = Math.abs((tstar - tstar_old) / tstar);
                } else {
                    testt = Math.abs(tstar - tstar_old);
                }

                double testq;
                if (qstar != 0.0) {
                    testq = Math.abs((qstar - qstar_old) / qstar);
                } else {
                    testq = Math.abs(qstar - qstar_old);
                }

                if (testu < 0.001 && testt < 0.001 && testq < 0.001) {
                    break;
                }
            }

            double cp_air = EvapUtilities.spec_ht_air(ts);
            double lv = EvapUtilities.latent(ts);
            EvapUtilities.DoubleContainer taudc = new EvapUtilities.DoubleContainer();
            EvapUtilities.DoubleContainer hsdc = new EvapUtilities.DoubleContainer();
            EvapUtilities.DoubleContainer hldc = new EvapUtilities.DoubleContainer();
            EvapUtilities.DoubleContainer evapdc = new EvapUtilities.DoubleContainer();
            EvapUtilities.fluxes(ustar, tstar, qstar, rhoa, cp_air, lv, taudc, hsdc, hldc, evapdc);
            double tau = taudc.d;
            this._hl = hldc.d;
            this._hs = hsdc.d;
            this._evap = evapdc.d;
            this._ustar = ustar;
            this._qstar = qstar;
            this._rstar = rstar;
            this._obukhovLen = l;
        }
    }
}
