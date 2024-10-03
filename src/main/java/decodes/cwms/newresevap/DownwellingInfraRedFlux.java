/* 
 * Copyright (c) 2018
 * United States Army Corps of Engineers - Hydrologic Engineering Center (USACE/HEC)
 * All Rights Reserved.  USACE PROPRIETARY/CONFIDENTIAL.
 * Source may not be released without written approval from HEC
 */
package decodes.cwms.newresevap;

import decodes.cwms.HecConstants;
import decodes.cwms.resevapcalc.CloudCover;
import decodes.cwms.resevapcalc.Const;

/**
 * Data and functions for computing downwelling IR flux. 

 * @author RESEVAP program by Steven F. Daly  (ERDC/CRREL)
 * conversion to Java by Richard Rachiele (RMA)
 */
public class DownwellingInfraRedFlux
{
    private static final double[][] COEF11 = new double[][]{{1.05, 0.6, 5.0, 25.0}, {4.1, 0.3, 4.0, 25.0}, {7.0, 1.5, 3.0, 30.0}};
    private static final double[][] COEF12 = new double[][]{{1.05, 0.6, 1.5, 25.0}, {4.1, 2.0, 1.7, 25.0}, {7.0, 1.5, 3.0, 30.0}};
    private static final double[][] COEF21 = new double[][]{{1.15, 0.45, 5.0, 25.0}, {4.1, 2.0, 1.7, 25.0}, {7.0, 1.5, 3.0, 30.0}};
    private static final double[][] COEF22 = new double[][]{{1.15, 0.6, 1.5, 25.0}, {4.4, 1.2, 3.0, 25.0}, {7.0, 1.5, 3.0, 30.0}};
    private static final double[][][][] COEF = new double[2][2][3][4];

    public Dnirflx() {
    }

    public static double dnirflx(int jday, double airTemp, double rh, double ematm, double lat, CloudCover[] cloudCover) {
        double doy = (double)jday;
        double ta = airTemp;
        double lcldbse = cloudCover[2].height;
        double mcldbse = cloudCover[1].height;
        double hcldbse = cloudCover[0].height;
        double lcld = cloudCover[2].fractionCloudCover;
        double mcld = cloudCover[1].fractionCloudCover;
        double hcld = cloudCover[0].fractionCloudCover;
        if (!RMAConst.isValidValue(lcld)) {
            lcld = cloudCover[2].getDefaultFractionCloudCover();
        }

        if (!RMAConst.isValidValue(mcld)) {
            mcld = cloudCover[1].getDefaultFractionCloudCover();
        }

        if (!RMAConst.isValidValue(hcld)) {
            hcld = cloudCover[0].getDefaultFractionCloudCover();
        }

        int ilat = 0;
        if (lat >= Math.abs(25.0)) {
            ilat = 1;
        }

        int isean = 1;
        if (lat > 0.0 && (int)doy > 330 || (int)doy < 65) {
            isean = 0;
        }

        if (lat < 0.0 && (int)doy > 150 || (int)doy < 250) {
            isean = 0;
        }

        double a;
        double b;
        double c;
        double d;
        double zlcld;
        if (!RMAConst.isValidValue(lcldbse) && lcld != 0.0) {
            a = COEF[isean][ilat][0][0];
            b = COEF[isean][ilat][0][1];
            c = COEF[isean][ilat][0][2];
            d = COEF[isean][ilat][0][3];
            zlcld = a - b * (1.0 - Math.abs(Math.cos(c * (lat - d))));
        } else {
            zlcld = lcldbse;
        }

        double zmcld;
        if (!RMAConst.isValidValue(mcldbse) && mcld != 0.0) {
            a = COEF[isean][ilat][1][0];
            b = COEF[isean][ilat][1][1];
            c = COEF[isean][ilat][1][2];
            d = COEF[isean][ilat][1][3];
            zmcld = a - b * (1.0 - Math.abs(Math.cos(c * (lat - d))));
        } else {
            zmcld = mcldbse;
        }

        double zhcld;
        if (!RMAConst.isValidValue(hcldbse) && hcld != 0.0) {
            a = COEF[isean][ilat][2][0];
            b = COEF[isean][ilat][2][1];
            c = COEF[isean][ilat][2][2];
            d = COEF[isean][ilat][2][3];
            zhcld = a - b * (1.0 - Math.abs(Math.cos(c * (lat - d))));
        } else {
            zhcld = hcldbse;
        }

        hcld = hcld * (1.0 - mcld) * (1.0 - lcld);
        mcld *= 1.0 - lcld;
        if (lcld == 0.0) {
            zlcld = 0.0;
        }

        if (mcld == 0.0) {
            zmcld = 0.0;
        }

        if (hcld == 0.0) {
            zhcld = 0.0;
        }

        double flxclr = ematm * 5.669E-8 * Math.pow(ta + 273.15, 4.0);
        double flxcld = lcld * (94.0 - 5.8 * zlcld) + mcld * (94.0 - 5.8 * zmcld) + hcld * (94.0 - 5.8 * zhcld);
        double flxir = flxclr + flxcld;
        return flxir;
    }

    public static double emisatm(double airTemp, double relH) {
        double rh = relH;
        double rv = 461.0;
        double eso = 6.13;
        double ta = airTemp + 273.15;
        double latent = (-0.00243 * ta + 3.166659) * 1000000.0;
        double ea = eso * Math.exp(latent / rv * (0.0036609921288669233 - 1.0 / ta)) * rh * 0.01;
        double ematm = 1.24 * Math.pow(ea / ta, 0.14285714285714285);
        return ematm;
    }

    static {
        COEF[0][0] = COEF11;
        COEF[0][1] = COEF12;
        COEF[1][0] = COEF21;
        COEF[1][1] = COEF22;
    }
    
}
