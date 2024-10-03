/* 
 * Copyright (c) 2018
 * United States Army Corps of Engineers - Hydrologic Engineering Center (USACE/HEC)
 * All Rights Reserved.  USACE PROPRIETARY/CONFIDENTIAL.
 * Source may not be released without written approval from HEC
 */
package decodes.cwms.newresevap;

//import hec.heclib.util.HecTime;
//import rma.util.RMAConst;

import decodes.cwms.HecConstants;
import decodes.cwms.resevapcalc.CloudCover;
import decodes.cwms.resevapcalc.Const;

import java.util.Calendar;
import java.util.Date;

/**
 *
 * @author RESEVAP program by Steven F. Daly (ERDC/CRREL)
 * conversion to Java by Richard Rachiele (RMA)
 */
public class SolarFlux
{
    private static final double[][] R1 = new double[][]{{0.12395, 0.15325, 0.15946, 0.27436}, {-0.34765, -0.3962, -0.42185, -0.43132}, {0.39478, 0.42095, 0.488, 0.2692}, {-0.14627, -0.142, -0.18492, -0.00447}};
    private static final double[][] R2 = new double[][]{{0.25674, 0.42111, 0.61394, 0.69143}, {-0.18077, -0.04002, -0.01469, -0.14419}, {-0.21961, -0.51833, -0.174, -0.051}, {0.25272, 0.4054, 0.14215, 0.06682}};
    private static final double[][] T1 = new double[][]{{0.76977, 0.69318, 0.68679, 0.55336}, {0.49407, 0.68227, 0.71012, 0.61511}, {-0.44647, -0.64289, -0.71463, -0.29816}, {0.11558, 0.1791, 0.22339, -0.06663}};
    private static final double[][] T2 = new double[][]{{0.63547, 0.43562, 0.23865, 0.15785}, {0.35229, 0.26094, 0.20143, 0.3241}, {0.08709, 0.36428, -0.01183, -0.14458}, {-0.22902, -0.38556, -0.07892, 0.01457}};
    private static final double[][] WTX = new double[][]{{0.675, 1.552, 1.429, 1.512}, {-3.432, -1.957, -1.207, -1.176}, {1.929, -1.762, -2.008, -2.16}, {0.842, 2.067, 0.853, 1.42}, {2.693, 0.448, 0.324, -0.032}, {-1.354, 0.932, 1.582, 1.422}};
    double _diffuse;
    double _direct;
    double _sdown;
    double ZEN;
    double SOLAR;

    public Solflx() {
    }

    public void solflx(HecTime currentTime, double gmtOffset, double longitude, double latitude, CloudCover[] cloudCover) {
        EvapUtilities.DoubleContainer solzen = new EvapUtilities.DoubleContainer();
        new EvapUtilities.DoubleContainer();
        EvapUtilities.DoubleContainer cosz = new EvapUtilities.DoubleContainer();
        int doy = currentTime.dayOfYear();
        int hr = currentTime.hour();
        double frad = 0.017453292;
        double dhr = (double)hr;
        if (longitude <= -180.0) {
            longitude += 360.0;
        }

        if (longitude > 180.0) {
            longitude -= 360.0;
        }

        zenith(doy, dhr, gmtOffset, latitude, longitude, solzen, cosz);
        double sdown;
        if (solzen.d >= 90.0) {
            sdown = 0.0;
            double direct = 0.0;
            double var20 = 0.0;
        } else {
            int jday = doy;
            this.insol(jday, cloudCover, cosz.d);
            sdown = this._sdown;
        }

        this.SOLAR = sdown;
        this.ZEN = solzen.d;
    }

    public void insol(int jday, CloudCover[] cloudCover, double cosz) {
        double[] rk = new double[3];
        double[] tk = new double[3];
        double[] tdk = new double[3];
        int[] icla = new int[3];
        double[] covera = new double[3];
        double sdowne = 6.283185308 * (double)((float)(jday - 2)) / 365.242;
        sdowne = 1.0001399 + 0.0167261 * Math.cos(sdowne);
        sdowne *= sdowne;
        double coszsq = cosz * cosz;
        double coszcube = coszsq * cosz;
        double sdown0 = 1369.2 * sdowne * cosz;
        covera[0] = cloudCover[0].fractionCloudCover;
        covera[1] = cloudCover[1].fractionCloudCover;
        covera[2] = cloudCover[2].fractionCloudCover;
        icla[0] = cloudCover[0].getCloudTypeFlag();
        icla[1] = cloudCover[1].getCloudTypeFlag();
        icla[2] = cloudCover[2].getCloudTypeFlag();
        if (!RMAConst.isValidValue(covera[2])) {
            covera[2] = cloudCover[2].getDefaultFractionCloudCover();
        }

        if (!RMAConst.isValidValue(covera[1])) {
            covera[1] = cloudCover[1].getDefaultFractionCloudCover();
        }

        if (!RMAConst.isValidValue(covera[0])) {
            covera[0] = cloudCover[0].getDefaultFractionCloudCover();
        }

        for(int i = 0; i < 3; ++i) {
            double wgt = 0.0;
            double rcld = 0.0;
            double tcld = 0.0;
            double rclr = 0.0;
            double tclr = 0.0;
            int jj = icla[i];
            double fr = covera[i];
            int ll = i;
            if (jj != 0) {
                int j = jj - 1;
                wgt = WTX[0][j] + WTX[1][j] * cosz + WTX[2][j] * fr + WTX[3][j] * cosz * fr + WTX[4][j] * coszsq + WTX[5][j] * fr * fr;
                wgt *= fr;
                if (fr < 0.05) {
                    wgt = 0.0;
                }

                if (fr > 0.95) {
                    wgt = 1.0;
                }

                if (wgt > 0.0) {
                    rcld = R2[0][j] + R2[1][j] * cosz + R2[2][j] * coszsq + R2[3][j] * coszcube;
                    tcld = T2[0][j] + T2[1][j] * cosz + T2[2][j] * coszsq + T2[3][j] * coszcube;
                }
            }

            if (wgt < 1.0) {
                rclr = R1[0][ll] + R1[1][ll] * cosz + R1[2][ll] * coszsq + R1[3][ll] * coszcube;
                tclr = T1[0][ll] + T1[1][ll] * cosz + T1[2][ll] * coszsq + T1[3][ll] * coszcube;
            }

            rk[i] = wgt * rcld + (1.0 - wgt) * rclr;
            tk[i] = wgt * tcld + (1.0 - wgt) * tclr;
            tdk[i] = tk[i] - rk[i];
            if (tdk[i] < 0.0) {
                tdk[i] = 0.0;
            }
        }

        double rg = 0.2;
        double d1 = 1.0 - rk[0] * rk[1];
        double d2 = 1.0 - rk[1] * rk[2];
        double d3 = 1.0 - rk[2] * rg;
        double sdown = d1 * d2 - rk[0] * rk[2] * tk[1] * tk[1];
        sdown = d3 * sdown - d1 * rk[1] * rg * tk[2] * tk[2];
        sdown -= rk[0] * rg * tk[1] * tk[1] * tk[2] * tk[2];
        sdown = tk[0] * tk[1] * tk[2] * sdown0 / sdown;
        double direct;
        double diffuse;
        if (sdown <= 0.0) {
            sdown = 0.0;
            direct = 0.0;
            diffuse = 0.0;
        } else {
            direct = tdk[0] * tdk[1] * tdk[2] * sdown0;
            diffuse = sdown - direct;
        }

        this._sdown = sdown;
        this._direct = direct;
        this._diffuse = diffuse;
    }

    public static void zenith(int jjday, double local_hr, double gmt_offset, double zlat, double zlong, EvapUtilities.DoubleContainer solzen, EvapUtilities.DoubleContainer coszd) {
        double gmt = local_hr + gmt_offset;
        double frad = 0.01745329252;
        double xpi = 180.0 * frad;
        double phi = 360.0 * (double)((float)(jjday - 1)) / 365.242;
        double phir = phi * frad;
        double sinp = Math.sin(phir);
        double cosp = Math.cos(phir);
        double sin2p = 2.0 * sinp * cosp;
        double cos2p = 2.0 * cosp * cosp - 1.0;
        double sig = 279.9348 + phi + 1.914827 * sinp - 0.079525 * cosp + 0.019938 * sin2p - 0.001639 * cos2p;
        double sigr = sig * frad;
        double sind = 0.39785 * Math.sin(sigr);
        double cosd = Math.sqrt(1.0 - sind * sind);
        double xm = 12.0 + 0.12357 * sinp - 0.004289 * cosp + 0.153809 * sin2p + 0.060783 * cos2p;
        double h = 15.0 * (gmt - xm) - zlong;
        double hr = h * frad;
        double zlatr = zlat * frad;
        double cosz = Math.sin(zlatr) * sind + Math.cos(zlatr) * cosd * Math.cos(hr);
        double sinz = Math.sqrt(1.0 - cosz * cosz);
        double saz;
        if (sinz == 0.0) {
            saz = xpi;
        } else if (sind - Math.sin(zlatr) * cosz <= 0.0) {
            double az = cosd * Math.sin(hr) / sinz;
            if (az < -1.0) {
                az = -1.0;
            }

            saz = Math.asin(az) + xpi;
        } else if (hr > 0.0) {
            saz = 2.0 * xpi - Math.asin(cosd * Math.sin(hr) / sinz);
        } else {
            saz = -Math.asin(cosd * Math.sin(hr) / sinz);
        }

        saz /= frad;
        if (saz < 0.0) {
            saz += 360.0;
        }

        solzen.d = Math.acos(cosz) / 0.017453292;
        coszd.d = cosz;
    }
}
