/* 
 * Copyright (c) 2018
 * United States Army Corps of Engineers - Hydrologic Engineering Center (USACE/HEC)
 * All Rights Reserved.  USACE PROPRIETARY/CONFIDENTIAL.
 * Source may not be released without written approval from HEC
 */
package decodes.cwms.resevapcalc;

import decodes.db.Constants;
import decodes.tsdb.CTimeSeries;
//import hec.heclib.util.HecTime;
//import hec.io.TimeSeriesContainer;
//import rma.util.RMAConst;

import java.util.Date;

/**
 *  Class hold meteorological time series data for Resevap program.
 */
public class EvapMetData
{
    // input met data timeseries
    CTimeSeries windspeedTsc;
    CTimeSeries airTempTsc;
    CTimeSeries relHumidityTsc;
    CTimeSeries dewPointTsc;
    CTimeSeries airPressureTsc;
    CTimeSeries fractionLowClouds;
    CTimeSeries altitudeLowClouds;
    CTimeSeries fractionMedClouds;
    CTimeSeries altitudeMedClouds;
    CTimeSeries fractionHighClouds;
    CTimeSeries altitudeHighClouds;
    
    // these variables hold last valid values
    //was originally RMA undefinedDouble = -FloatMax
    double wsTemp_old = Constants.undefinedDouble;
    double windSpeed_old = Constants.undefinedDouble;
    double airTemp_old = Constants.undefinedDouble;
    double relHumidity_old = Constants.undefinedDouble;
    double airPressure_old = Constants.undefinedDouble;
    
    // test variables hold the met values for the current time
    double wsTemp_current = Constants.undefinedDouble;
    double windSpeed_current = Constants.undefinedDouble;
    double airTemp_current = Constants.undefinedDouble;
    double relHumidity_current = Constants.undefinedDouble;
    double airPressure_current = Constants.undefinedDouble;
          
    public double getWindSpeed( Date hecTime ) throws ResEvapException {
        return getMetValue(windspeedTsc, hecTime);
    }
    
    public double getAirTemp( Date hecTime ) throws ResEvapException {
        return getMetValue(airTempTsc, hecTime);
    }
    public double getRelHumidity( Date hecTime ) throws ResEvapException {
        return getMetValue(relHumidityTsc, hecTime);
    }
    
    public double getDewPoint( Date hecTime ) throws ResEvapException {
        return getMetValue(dewPointTsc, hecTime);
    }
    
    public double getAirPressure( Date hecTime ) throws ResEvapException {
        return getMetValue(airPressureTsc, hecTime);
    }
    
    /**
     * Get CloudCover array for current time 
     * These hold the fractional cloud cover and 
     * base height for the Low, Mid and High cloud
     * divisions.
     * 
     * @param hecTime
     * @return 
     */
    public CloudCover[] getCloudCover( Date hecTime ) throws ResEvapException {
        CloudCover[] cloudCover = new CloudCover[3];
        double fractionCC = getMetValue(fractionLowClouds, hecTime);
        double altitude = getMetValue(altitudeLowClouds, hecTime);
        
        cloudCover[2] = new CloudCover( fractionCC,
                altitude, CloudCover.CloudHeightType.height_low);

        fractionCC = getMetValue(fractionMedClouds, hecTime);
        altitude = getMetValue(altitudeMedClouds, hecTime);
        
        cloudCover[1] = new CloudCover( fractionCC,
                    altitude, CloudCover.CloudHeightType.height_med);

        fractionCC = getMetValue(fractionHighClouds, hecTime);
        altitude = getMetValue(altitudeHighClouds, hecTime);
        
        cloudCover[0] = new CloudCover( fractionCC,
                altitude, CloudCover.CloudHeightType.height_high);
        
        return cloudCover;
    }

    private double getMetValue( CTimeSeries tsc, Date hecTime ) throws ResEvapException {
            int idx = tsc.findNextIdx(hecTime);
            double value;
            try{
                value = tsc.sampleAt(idx).getDoubleValue();
            }
            catch (Exception ex){
                throw new ResEvapException("failed to load met value from timeseries", ex);
            }
            return value;
    }
    
    public void setAirTempTs( CTimeSeries tsc )
    {
        airTempTsc = tsc;
    }
    
    public void setAirPressureTs( CTimeSeries tsc )
    {
        airPressureTsc = tsc;
    }

    public void setRelHumidityTs( CTimeSeries tsc )
    {
        relHumidityTsc = tsc;
    }

    public void setWindSpeedTs( CTimeSeries tsc )
    {
        windspeedTsc = tsc;
    }
    
    public void setHighCloudTs( CTimeSeries tscFrac, CTimeSeries tscHeight )
    {
        fractionHighClouds = tscFrac;
        altitudeHighClouds = tscHeight;
    }
    
    public void setMedCloudTs( CTimeSeries tscFrac, CTimeSeries tscHeight )
    {
        fractionMedClouds = tscFrac;
        altitudeMedClouds = tscHeight;
    }
    
    public void setLowCloudTs( CTimeSeries tscFrac, CTimeSeries tscHeight )
    {
        fractionLowClouds = tscFrac;
        altitudeLowClouds = tscHeight;
    }

}
