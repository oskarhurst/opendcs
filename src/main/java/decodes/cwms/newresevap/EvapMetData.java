/* 
 * Copyright (c) 2018
 * United States Army Corps of Engineers - Hydrologic Engineering Center (USACE/HEC)
 * All Rights Reserved.  USACE PROPRIETARY/CONFIDENTIAL.
 * Source may not be released without written approval from HEC
 */
package decodes.cwms.newresevap;

import decodes.cwms.resevapcalc.CloudCover;
import decodes.cwms.resevapcalc.ResEvapException;
import decodes.db.Constants;
import decodes.tsdb.CTimeSeries;

import java.util.Date;

/**
 *  Class hold meteorological time series data for Resevap program.
 */
public class EvapMetData
{
    TimeSeriesContainer _windspeedTsc;
    TimeSeriesContainer _airTempTsc;
    TimeSeriesContainer _relHumidityTsc;
    TimeSeriesContainer _dewPointTsc;
    TimeSeriesContainer _airPressureTsc;
    TimeSeriesContainer _fractionLowClouds;
    TimeSeriesContainer _altitudeLowClouds;
    TimeSeriesContainer _fractionMedClouds;
    TimeSeriesContainer _altitudeMedClouds;
    TimeSeriesContainer _fractionHighClouds;
    TimeSeriesContainer _altitudeHighClouds;
    double _wsTemp_old = Double.NEGATIVE_INFINITY;
    double _windSpeed_old = Double.NEGATIVE_INFINITY;
    double _airTemp_old = Double.NEGATIVE_INFINITY;
    double _relHumidity_old = Double.NEGATIVE_INFINITY;
    double _airPressure_old = Double.NEGATIVE_INFINITY;
    double _wsTemp_current = Double.NEGATIVE_INFINITY;
    double _windSpeed_current = Double.NEGATIVE_INFINITY;
    double _airTemp_current = Double.NEGATIVE_INFINITY;
    double _relHumidity_current = Double.NEGATIVE_INFINITY;
    double _airPressure_current = Double.NEGATIVE_INFINITY;

    public EvapMetData() {
    }

    public double getWindSpeed(HecTime hecTime) {
        return this.getMetValue(this._windspeedTsc, hecTime);
    }

    public double getAirTemp(HecTime hecTime) {
        return this.getMetValue(this._airTempTsc, hecTime);
    }

    public double getRelHumidity(HecTime hecTime) {
        return this.getMetValue(this._relHumidityTsc, hecTime);
    }

    public double getDewPoint(HecTime hecTime) {
        return this.getMetValue(this._dewPointTsc, hecTime);
    }

    public double getAirPressure(HecTime hecTime) {
        return this.getMetValue(this._airPressureTsc, hecTime);
    }

    public CloudCover[] getCloudCover(HecTime hecTime) {
        CloudCover[] cloudCover = new CloudCover[3];
        double fractionCC = this.getMetValue(this._fractionLowClouds, hecTime);
        double altitude = this.getMetValue(this._altitudeLowClouds, hecTime);
        cloudCover[2] = new CloudCover(fractionCC, altitude, CloudHeightType.height_low);
        fractionCC = this.getMetValue(this._fractionMedClouds, hecTime);
        altitude = this.getMetValue(this._altitudeMedClouds, hecTime);
        cloudCover[1] = new CloudCover(fractionCC, altitude, CloudHeightType.height_med);
        fractionCC = this.getMetValue(this._fractionHighClouds, hecTime);
        altitude = this.getMetValue(this._altitudeHighClouds, hecTime);
        cloudCover[0] = new CloudCover(fractionCC, altitude, CloudHeightType.height_high);
        return cloudCover;
    }

    private double getMetValue(TimeSeriesContainer tsc, HecTime hecTime) {
        if (tsc.interval > 0) {
            return tsc.getValue(hecTime);
        } else {
            long findTime = (long)hecTime.value();
            long sTimeValue = (long)tsc.startTime;
            long eTimeValue = (long)tsc.endTime;
            if (findTime >= sTimeValue && findTime <= eTimeValue) {
                for(int i = 0; i < tsc.times.length; ++i) {
                    if ((long)tsc.times[i] == findTime) {
                        return tsc.values[i];
                    }

                    if ((long)tsc.times[i] > findTime) {
                        return -3.4028234663852886E38;
                    }
                }

                return -3.4028234663852886E38;
            } else {
                return -3.4028234663852886E38;
            }
        }
    }

    public void setAirTempTs(TimeSeriesContainer tsc) {
        this._airTempTsc = tsc;
    }

    public void setAirPressureTs(TimeSeriesContainer tsc) {
        this._airPressureTsc = tsc;
    }

    public void setRelHumidityTs(TimeSeriesContainer tsc) {
        this._relHumidityTsc = tsc;
    }

    public void setWindSpeedTs(TimeSeriesContainer tsc) {
        this._windspeedTsc = tsc;
    }

    public void setHighCloudTs(TimeSeriesContainer tscFrac, TimeSeriesContainer tscHeight) {
        this._fractionHighClouds = tscFrac;
        this._altitudeHighClouds = tscHeight;
    }

    public void setMedCloudTs(TimeSeriesContainer tscFrac, TimeSeriesContainer tscHeight) {
        this._fractionMedClouds = tscFrac;
        this._altitudeMedClouds = tscHeight;
    }

    public void setLowCloudTs(TimeSeriesContainer tscFrac, TimeSeriesContainer tscHeight) {
        this._fractionLowClouds = tscFrac;
        this._altitudeLowClouds = tscHeight;
    }
}
