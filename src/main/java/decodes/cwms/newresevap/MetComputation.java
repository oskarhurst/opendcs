/* 
 * Copyright (c) 2018
 * United States Army Corps of Engineers - Hydrologic Engineering Center (USACE/HEC)
 * All Rights Reserved.  USACE PROPRIETARY/CONFIDENTIAL.
 * Source may not be released without written approval from HEC
 */
package decodes.cwms.newresevap;

import decodes.cwms.HecConstants;
import decodes.cwms.resevapcalc.CloudCover;
import decodes.cwms.resevapcalc.DownwellingInfraRedFlux;
import decodes.cwms.resevapcalc.EvapMetData;

import java.util.Calendar;
import java.util.Date;
import java.util.logging.Level;
import java.util.logging.Logger;
//import rma.util.RMAConst;

/**
 * Compute surface heat exchange and evaporation rate from
 * meteorological values.
 * 
 * @author richard
 */
public class MetComputation
{
    private static final Logger LOGGER = Logger.getLogger(MetComputation.class.getName());
    EvapMetData _metData;
    EvapWater _evapWater = new EvapWater();
    double _solar;
    double _flxir;
    double _flxir_out;
    double _hs;
    double _hl;
    boolean _metFailed;
    boolean _metDefined = false;

    public MetComputation() {
    }

    public void setMetData(EvapMetData metData) {
        this._metData = metData;
    }

    public void setEvapWater(EvapWater evapWater) {
        this._evapWater = evapWater;
    }

    public void computeMetAndEvap(HecTime currentTime, double surfaceTemp, ReservoirLocationInfo resLocationInfo) {
        double airPressure = this._metData.getAirPressure(currentTime);
        if (airPressure < 0.0) {
            airPressure = -901.0;
        }

        double windSpeed = this._metData.getWindSpeed(currentTime);
        double relHumidity = this._metData.getRelHumidity(currentTime);
        double airTemp = this._metData.getAirTemp(currentTime);
        CloudCover[] cloudCover = this._metData.getCloudCover(currentTime);
        boolean missingMetData = false;
        this._metFailed = false;
        String msg;
        if (!RMAConst.isValidValue(windSpeed)) {
            msg = "Wind speed missing or bad " + currentTime;
            LOGGER.log(Level.SEVERE, msg);
            missingMetData = true;
        }

        if (!RMAConst.isValidValue(airTemp)) {
            msg = "Air temp missing or bad " + currentTime;
            LOGGER.log(Level.SEVERE, msg);
            missingMetData = true;
        }

        if (!RMAConst.isValidValue(relHumidity)) {
            msg = "RH missing or bad " + currentTime;
            LOGGER.log(Level.SEVERE, msg);
            missingMetData = true;
        }

        if (!RMAConst.isValidValue(airPressure)) {
            msg = "PRESSURE missing or bad " + currentTime;
            LOGGER.log(Level.SEVERE, msg);
            missingMetData = true;
        }

        if (!missingMetData) {
            this._metDefined = true;
        }

        if (!this._metDefined) {
            this._metFailed = true;
        } else {
            boolean cloudFailed = false;

            for(int ic = 0; ic < 3; ++ic) {
                CloudCover cld = cloudCover[ic];
                String msg;
                if (!RMAConst.isValidValue(cld.fractionCloudCover)) {
                    int ic1 = ic + 1;
                    msg = " Cover (" + ic1 + ") missing or bad" + currentTime;
                    LOGGER.log(Level.SEVERE, msg);
                }

                if (!RMAConst.isValidValue(cld.height)) {
                    String heightStr = cld.getTypeName();
                    msg = heightStr + " missing or bad " + currentTime;
                    LOGGER.log(Level.SEVERE, msg);
                }
            }

            if (!RMAConst.isValidValue(windSpeed)) {
                windSpeed = this._metData._windSpeed_old;
            } else {
                this._metData._windSpeed_old = windSpeed;
            }

            if (!RMAConst.isValidValue(airTemp)) {
                airTemp = this._metData._airTemp_old;
            } else {
                this._metData._airTemp_old = airTemp;
            }

            if (!RMAConst.isValidValue(relHumidity)) {
                relHumidity = this._metData._relHumidity_old;
            } else {
                this._metData._relHumidity_old = relHumidity;
            }

            if (!RMAConst.isValidValue(airPressure)) {
                airPressure = this._metData._airPressure_old;
            } else {
                this._metData._airPressure_old = airPressure;
            }

            this._metData._windSpeed_current = windSpeed;
            this._metData._airTemp_current = airTemp;
            this._metData._relHumidity_current = relHumidity;
            this._metData._airPressure_current = airPressure;
            double longitude = resLocationInfo.lon;
            double latitude = resLocationInfo.lat;
            double gmtOffset = resLocationInfo.gmtOffset;
            this.computeMet(currentTime, gmtOffset, surfaceTemp, windSpeed, airTemp, airPressure, relHumidity, cloudCover, latitude, longitude);
            if (windSpeed < 0.1) {
                windSpeed = 0.1;
                this._metData._windSpeed_current = windSpeed;
            }

            if (surfaceTemp < 0.0) {
                surfaceTemp = 0.0;
            }

            this.computeEvap(resLocationInfo, surfaceTemp, windSpeed, airTemp, airPressure, relHumidity);
        }
    }

    public void computeMet(HecTime currentTime, double gmtOffset, double surfaceTemp, double windSpeed, double airTemp, double airPressure, double relHumidity, CloudCover[] cloudCover, double lat, double lon) {
        int jday = currentTime.dayOfYear();
        CloudCover[] var20 = cloudCover;
        int var21 = cloudCover.length;

        for(int var22 = 0; var22 < var21; ++var22) {
            CloudCover cloud = var20[var22];
            if (RMAConst.isValidValue(cloud.height)) {
                cloud.height /= 1000.0;
            }
        }

        Solflx solarFlx = new Solflx();
        solarFlx.solflx(currentTime, gmtOffset, lon, lat, cloudCover);
        this._solar = solarFlx.SOLAR;
        double ws_temp;
        double flxir_out;
        if (RMAConst.isValidValue(surfaceTemp)) {
            ws_temp = Dnirflx.emisatm(airTemp, relHumidity);
            flxir_out = Dnirflx.dnirflx(jday, airTemp, relHumidity, ws_temp, lat, cloudCover);
            this._flxir = flxir_out;
        }

        ws_temp = surfaceTemp + 273.15;
        flxir_out = -5.5565999999999996E-8 * Math.pow(ws_temp, 4.0);
        this._flxir_out = flxir_out;
    }

    public void computeEvap(ReservoirLocationInfo resLocationInfo, double surfaceTemp, double windSpeed, double airTemp, double airPressure, double relHumidity) {
        double rt = resLocationInfo.rt;
        double ru = resLocationInfo.ru;
        double rq = resLocationInfo.rq;
        if (this._evapWater == null) {
            this._evapWater = new EvapWater();
        }

        this._evapWater.evap_water(surfaceTemp, airTemp, relHumidity, windSpeed, airPressure, rt, ru, rq);
    }
       
}
