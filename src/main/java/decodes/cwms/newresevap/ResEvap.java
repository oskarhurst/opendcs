/* 
 * Copyright (c) 2018
 * United States Army Corps of Engineers - Hydrologic Engineering Center (USACE/HEC)
 * All Rights Reserved.  USACE PROPRIETARY/CONFIDENTIAL.
 * Source may not be released without written approval from HEC
 */
package decodes.cwms.newresevap;

import decodes.cwms.resevapcalc.EvapMetData;
import decodes.cwms.resevapcalc.EvapReservoir;
import decodes.tsdb.CTimeSeries;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;

//import rma.lang.RmaMath;
//import rma.util.RMAConst;

/**
 * Program to estimate evaporation from reservoirs
 * Program estimates 1-D water temperature profile in reservoir  
 * Calculates evaporation
 * 
 * @author RESEVAP program by Steven F. Daly (ERDC/CRREL)
 * Version 1.1  31 August 2009 
 * VERSION 2.0 7 July 2010
 * conversion to Java by Richard Rachiele (RMA)
 */
public class ResEvap
{
    private static final Logger LOGGER = Logger.getLogger(ResEvap.class.getName());
    public static final double EMITTANCE_H20 = 0.98;
    public static final double PENFRAC = 0.4;
    public static final double ALBEDO = 0.08;
    public static final double THETA = 1.0;
    public static final double SIGMA = 5.67E-8;
    public static final double ETA_CONVECTIVE = 0.5;
    public static final double ETA_STIRRING = 0.4;
    public static final double WIND_CRITIC = 1.0;
    NumberFormat nf8_3 = new NumberFormat("%8.3f");
    NumberFormat nf8_2 = new NumberFormat("%8.2f");
    NumberFormat nf9_3 = new NumberFormat("%9.3f");
    NumberFormat nf10_3 = new NumberFormat("%10.3f");
    NumberFormat nfe11_3 = new NumberFormat("%11.3E");
    TimeSeriesContainer _solarRadTsc;
    TimeSeriesContainer _IR_DownTsc;
    TimeSeriesContainer _IR_OutTsc;
    TimeSeriesContainer _sensibleHeatTsc;
    TimeSeriesContainer _latentHeatTsc;
    TimeSeriesContainer _evapRateHourlyTsc;
    TimeSeriesContainer _evapDailyTsc;
    TimeSeriesContainer _surfaceTempTsc;
    NavigableMap<Integer, Integer> _timeMap;
    TimeSeriesContainer[] _inputTimeSeries;
    String _versionName;
    double[][] _wtempProfiles;
    public EvapReservoir _reservoir;
    public EvapMetData _metData;
    private File _workDir;

    public ResEvap() {
    }

    public ResEvap(EvapReservoir reservoir, EvapMetData metData) {
        this.setReservoir(reservoir);
        this.setMetData(metData);
    }

    public boolean setOutputDirectory(File workDir) {
        this._workDir = workDir;
        return true;
    }

    public boolean setReservoir(EvapReservoir reservoir) {
        this._reservoir = reservoir;
        if (!reservoir.initRes()) {
            return false;
        } else {
            return reservoir.resSetup();
        }
    }

    public boolean setMetData(EvapMetData metData) {
        this._metData = metData;
        return true;
    }

    public boolean compute(String startDateTime, String endDateTime, double gmtOffset, String versionName) throws ResEvapException {
        if (this._reservoir == null) {
            throw new ResEvapException("ResEvap.compute: No reservoir has been set");
        } else if (this._metData == null) {
            throw new ResEvapException("ResEvap.compute: No Meteorological data has been set");
        } else {
            File outfil = null;
            File metoutfil = null;
            File toutfil = null;
            File xoutfil = null;
            if (this._workDir != null) {
                outfil = new File(this._workDir.getAbsolutePath() + "/wtout_java.dat");
                metoutfil = new File(this._workDir.getAbsolutePath() + "/testtout_java.dat");
                toutfil = new File(this._workDir.getAbsolutePath() + "/xout2_java.dat");
                xoutfil = new File(this._workDir.getAbsolutePath() + "/xout_java.dat");
            }

            this._versionName = versionName;
            this.initializeOutputTsc(this._reservoir._name, versionName, startDateTime, endDateTime);
            HecTime hecStartTime = new HecTime(startDateTime, 1);
            HecTime hecEndTime = new HecTime(endDateTime, 1);
            int intervalMinutes = 60;
            int nper = HecTime.nopers(intervalMinutes, hecStartTime.julian(), hecStartTime.minutesSinceMidnight(), hecEndTime.julian(), hecEndTime.minutesSinceMidnight());
            double deltT = 3600.0;
            MetComputation metComputation = new MetComputation();
            metComputation.setMetData(this._metData);
            ResWtCompute resWtCompute = new ResWtCompute(this._reservoir);
            BufferedWriter tout = null;
            if (toutfil != null) {
                try {
                    tout = new BufferedWriter(new FileWriter(toutfil));
                    resWtCompute.setOutfile(tout);
                } catch (IOException var46) {
                    IOException ex = var46;
                    LOGGER.log(Level.FINE, "Unable to read " + toutfil.getAbsolutePath(), ex);
                }
            }

            BufferedWriter xout = null;
            if (xoutfil != null) {
                try {
                    xout = new BufferedWriter(new FileWriter(xoutfil));
                    this._reservoir.setDebugFile(xout);
                } catch (IOException var45) {
                    IOException ex = var45;
                    LOGGER.log(Level.FINE, "Unable to read " + xoutfil.getAbsolutePath(), ex);
                }
            }

            BufferedWriter out = null;
            if (outfil != null) {
                try {
                    out = new BufferedWriter(new FileWriter(outfil));
                } catch (IOException var44) {
                    IOException ex = var44;
                    LOGGER.log(Level.FINE, "Unable to read " + outfil.getAbsolutePath(), ex);
                }
            }

            BufferedWriter metout = null;
            IOException ex;
            if (outfil != null) {
                try {
                    metout = new BufferedWriter(new FileWriter(metoutfil));
                    if (metout != null) {
                        String heading = "    Date    JD  GMT     U       T      RH       P       Ts      K       u*        R*       L          Hs        HL        Qs        IR       IR_out     Evap";
                        metout.write(heading);
                        metout.newLine();
                        heading = "                       m/s    deg C     %      mb      deg C           m/s                m      ********** W/m**2 ***********************     mm/d";
                        metout.write(heading);
                        metout.newLine();
                    }
                } catch (IOException var43) {
                    ex = var43;
                    LOGGER.log(Level.FINE, "Unable to read " + metoutfil.getAbsolutePath(), ex);
                }
            }

            try {
                if (xout != null) {
                    this._reservoir.resSetup();
                }

                int resj = this._reservoir.getResj();
                this._reservoir._resj_old = resj;
                HecTime currentTime = new HecTime(hecStartTime);
                ReservoirLocationInfo resLocationInfo = this._reservoir.getReservoirLocationInfo();
                resLocationInfo.gmtOffset = gmtOffset;
                boolean useElevTS = false;
                double wselCurrent = this._reservoir.getElevation();

                for(int jhour = 0; jhour <= nper; ++jhour) {
                    double surfaceTemp;
                    if (useElevTS) {
                        surfaceTemp = this._reservoir.getCurrentElevation(currentTime);
                        if (RMAConst.isValidValue(surfaceTemp)) {
                            wselCurrent = surfaceTemp;
                            if (xout != null) {
                                xout.write(currentTime.date(4) + " " + currentTime.getTime(false) + "  wselCurrent  " + (float)wselCurrent);
                            }

                            this._reservoir.setElevationMeters(wselCurrent);
                            this._reservoir.resSetup(true);
                        }
                    }

                    if (jhour > 0) {
                        currentTime.addMinutes(intervalMinutes);
                    }

                    surfaceTemp = this._reservoir._wt[resj];
                    metComputation.computeMetAndEvap(currentTime, surfaceTemp, resLocationInfo);
                    boolean noProblem = true;
                    if (!metComputation._metFailed) {
                        noProblem = resWtCompute.computeReservoirTemp(currentTime, metComputation, deltT);
                    }

                    if (noProblem) {
                        int i;
                        if (out != null) {
                            resj = this._reservoir.getResj();
                            String dateTimeStr = currentTime.date(4) + "            " + currentTime.getTime(false);

                            try {
                                out.write(dateTimeStr);

                                for(i = resj; i >= 0; --i) {
                                    String strval = this.nf9_3.form(this._reservoir._wt[i]);
                                    out.write(strval);
                                }

                                out.newLine();
                            } catch (IOException var47) {
                                IOException ex = var47;
                                LOGGER.log(Level.FINE, "Unable to read " + outfil.getAbsolutePath(), ex);
                            }

                            this.outputMetComputation(currentTime, this._metData, metComputation, surfaceTemp, metout);
                        }

                        int idx = (Integer)this._timeMap.get(currentTime.value());
                        if (idx >= 0 && idx < this._solarRadTsc.times.length) {
                            this._solarRadTsc.values[idx] = metComputation._solar;
                            this._IR_DownTsc.values[idx] = metComputation._flxir;
                            this._IR_OutTsc.values[idx] = metComputation._flxir_out;
                            this._latentHeatTsc.values[idx] = metComputation._evapWater._hl;
                            this._sensibleHeatTsc.values[idx] = metComputation._evapWater._hs;
                            this._surfaceTempTsc.values[idx] = surfaceTemp;
                            this._evapRateHourlyTsc.values[idx] = metComputation._evapWater._evap / 24.0;
                            i = this._reservoir.getResj() + 1;
                            this._wtempProfiles[idx] = new double[1000];

                            for(int ilyr = 0; ilyr < i; ++ilyr) {
                                this._wtempProfiles[idx][ilyr] = this._reservoir._wt[ilyr];
                            }
                        }
                    }
                }

                try {
                    if (out != null) {
                        out.close();
                    }

                    if (metout != null) {
                        metout.close();
                    }

                    if (tout != null) {
                        tout.close();
                    }

                    if (xout != null) {
                        xout.close();
                    }
                } catch (IOException var42) {
                    LOGGER.log(Level.SEVERE, "IOException occurred while closing files", var42);
                }

                return true;
            } catch (RuntimeException | IOException var48) {
                ex = var48;

                try {
                    if (out != null) {
                        out.close();
                    }

                    if (metout != null) {
                        metout.close();
                    }

                    if (tout != null) {
                        tout.close();
                    }

                    if (xout != null) {
                        xout.close();
                    }
                } catch (IOException var41) {
                    IOException ioe = var41;
                    LOGGER.log(Level.SEVERE, "IOException occurred while closing files", ioe);
                }

                LOGGER.log(Level.SEVERE, "Error within computation", ex);
                throw new ResEvapException(ex);
            }
        }
    }

    private void outputMetComputation(HecTime currentTime, EvapMetData metData, MetComputation metComputation, double surfaceTemp, BufferedWriter metout) {
        try {
            String dateTimeStr = currentTime.date(4) + " " + currentTime.dayOfYear() + " " + currentTime.getTime(false);
            metout.write(dateTimeStr);
            double val = metData._windSpeed_current;
            if (!RMAConst.isValidValue(val)) {
                val = -901.0;
            }

            String strval = this.nf8_2.form(val);
            metout.write(strval);
            val = metData._airTemp_current;
            if (!RMAConst.isValidValue(val)) {
                val = -901.0;
            }

            strval = this.nf8_2.form(val);
            metout.write(strval);
            val = metData._relHumidity_current;
            if (!RMAConst.isValidValue(val)) {
                val = -901.0;
            }

            strval = this.nf8_2.form(val);
            metout.write(strval);
            val = metData._airPressure_current;
            if (!RMAConst.isValidValue(val)) {
                val = -901.0;
            }

            strval = this.nf8_2.form(val);
            metout.write(strval);
            val = surfaceTemp;
            if (!RMAConst.isValidValue(val)) {
                val = -901.0;
            }

            strval = this.nf9_3.form(val);
            metout.write(strval);
            metout.write("  -901.");
            val = metComputation._evapWater._ustar;
            if (!RMAConst.isValidValue(val)) {
                val = -901.0;
            }

            strval = this.nf9_3.form(val);
            metout.write(strval);
            val = metComputation._evapWater._rstar;
            if (!RMAConst.isValidValue(val)) {
                val = -901.0;
            }

            strval = this.nf9_3.form(val);
            metout.write(strval);
            val = metComputation._evapWater._obukhovLen;
            if (!RMAConst.isValidValue(val)) {
                val = -901.0;
            }

            if (val == 0.0) {
                val = -901.0;
            }

            float fval = (float)val;
            strval = this.nfe11_3.form((double)fval);
            metout.write(strval);
            val = metComputation._evapWater._hs;
            if (!RMAConst.isValidValue(val)) {
                val = -901.0;
            }

            strval = this.nf10_3.form(-val);
            metout.write(strval);
            val = metComputation._evapWater._hl;
            if (!RMAConst.isValidValue(val)) {
                val = -901.0;
            }

            strval = this.nf10_3.form(-val);
            metout.write(strval);
            val = metComputation._solar;
            if (!RMAConst.isValidValue(val)) {
                val = -901.0;
            }

            strval = this.nf10_3.form(val);
            metout.write(strval);
            val = metComputation._flxir;
            if (!RMAConst.isValidValue(val)) {
                val = -901.0;
            }

            strval = this.nf10_3.form(val);
            metout.write(strval);
            val = metComputation._flxir_out;
            if (!RMAConst.isValidValue(val)) {
                val = -901.0;
            }

            strval = this.nf10_3.form(val);
            metout.write(strval);
            val = metComputation._evapWater._evap;
            if (!RMAConst.isValidValue(val)) {
                val = -901.0;
            }

            strval = this.nf9_3.form(val);
            metout.write(strval);
            metout.newLine();
        } catch (IOException var12) {
            IOException ioe = var12;
            LOGGER.log(Level.FINE, "Unable to write output to ", ioe);
        }

    }

    protected boolean initializeOutputTsc(String locationName, String versionName, String startDateTime, String endDateTime) {
        try {
            int intervalMinutes = 60;
            int offsetMinutes = 0;
            TimeSeriesMath tsMath = (TimeSeriesMath)TimeSeriesMath.generateRegularIntervalTimeSeries(startDateTime, endDateTime, intervalMinutes, offsetMinutes, -3.4028234663852886E38);
            TimeSeriesContainer hourlyTsc = tsMath.getContainer();
            hourlyTsc.location = locationName;
            hourlyTsc.type = "Inst";
            hourlyTsc.version = versionName;
            DSSPathString dsspath = new DSSPathString("", locationName, "", "", "1HOUR", versionName);
            this._solarRadTsc = (TimeSeriesContainer)hourlyTsc.clone();
            this._solarRadTsc.parameter = "Irrad-Flux-Solar";
            this._solarRadTsc.units = "W/m2";
            dsspath.setCPart(this._solarRadTsc.parameter);
            this._solarRadTsc.fullName = dsspath.toString();
            this._IR_DownTsc = (TimeSeriesContainer)hourlyTsc.clone();
            this._IR_DownTsc.parameter = "Irrad-Flux-IR";
            this._IR_DownTsc.units = "W/m2";
            dsspath.setCPart(this._IR_DownTsc.parameter);
            this._IR_DownTsc.fullName = dsspath.toString();
            this._IR_OutTsc = (TimeSeriesContainer)hourlyTsc.clone();
            this._IR_OutTsc.parameter = "Irrad-Flux-Out";
            this._IR_OutTsc.units = "W/m2";
            dsspath.setCPart(this._IR_OutTsc.parameter);
            this._IR_OutTsc.fullName = dsspath.toString();
            this._surfaceTempTsc = (TimeSeriesContainer)hourlyTsc.clone();
            this._surfaceTempTsc.parameter = "Temp-Water-Surface";
            this._surfaceTempTsc.units = "C";
            dsspath.setCPart(this._surfaceTempTsc.parameter);
            this._surfaceTempTsc.fullName = dsspath.toString();
            this._sensibleHeatTsc = (TimeSeriesContainer)hourlyTsc.clone();
            this._sensibleHeatTsc.parameter = "Irrad-Heat-Sensible";
            this._sensibleHeatTsc.units = "W/m2";
            dsspath.setCPart(this._sensibleHeatTsc.parameter);
            this._sensibleHeatTsc.fullName = dsspath.toString();
            this._latentHeatTsc = (TimeSeriesContainer)hourlyTsc.clone();
            this._latentHeatTsc.parameter = "Irrad-Heat-Latent";
            this._latentHeatTsc.units = "W/m2";
            dsspath.setCPart(this._latentHeatTsc.parameter);
            this._latentHeatTsc.fullName = dsspath.toString();
            this._evapRateHourlyTsc = (TimeSeriesContainer)hourlyTsc.clone();
            this._evapRateHourlyTsc.parameter = "EvapRate";
            this._evapRateHourlyTsc.units = "mm/hr";
            dsspath.setCPart(this._evapRateHourlyTsc.parameter);
            this._evapRateHourlyTsc.fullName = dsspath.toString();
            this._timeMap = new TreeMap();
            HecTime hectim = new HecTime(hourlyTsc.times[0], 1);

            int nprofiles;
            for(nprofiles = 0; nprofiles < hourlyTsc.times.length; ++nprofiles) {
                hectim.set(hourlyTsc.times[nprofiles]);
                this._timeMap.put(hourlyTsc.times[nprofiles], nprofiles);
            }

            nprofiles = hourlyTsc.times.length;
            this._wtempProfiles = new double[nprofiles][];
            return true;
        } catch (HecMathException var12) {
            HecMathException hme = var12;
            LOGGER.log(Level.SEVERE, "IOException occurred while closing files", hme);
            return false;
        }
    }

    public List<TimeSeriesContainer> getComputedMetTimeSeries() {
        List<TimeSeriesContainer> computedTsList = new ArrayList();
        computedTsList.add((TimeSeriesContainer)this._surfaceTempTsc.clone());
        computedTsList.add((TimeSeriesContainer)this._sensibleHeatTsc.clone());
        computedTsList.add((TimeSeriesContainer)this._latentHeatTsc.clone());
        computedTsList.add((TimeSeriesContainer)this._solarRadTsc.clone());
        computedTsList.add((TimeSeriesContainer)this._IR_DownTsc.clone());
        computedTsList.add((TimeSeriesContainer)this._IR_OutTsc.clone());
        computedTsList.add((TimeSeriesContainer)this._evapRateHourlyTsc.clone());
        return computedTsList;
    }

    public TimeSeriesContainer getHourlyEvapRateTimeSeries() {
        return this._evapRateHourlyTsc;
    }

    public TimeSeriesContainer getHourlyEvapTimeSeries() {
        TimeSeriesMath math = null;

        try {
            TimeSeriesMath tsMath = new TimeSeriesMath(this._evapRateHourlyTsc);
            math = (TimeSeriesMath)tsMath.transformTimeSeries("1HOUR", "", "AVE", false);
        } catch (HecMathException var3) {
            HecMathException ex = var3;
            LOGGER.log(Level.SEVERE, "HecMathException occurred while converting hourly evap rate to hourly total evap", ex);
        }

        TimeSeriesContainer output = null;
        if (math != null) {
            output = math.getContainer();
            output.units = "mm";
            output.parameter = "Evap";
            output.type = DssDataType.PER_CUM.toString();
        }

        return output;
    }

    public TimeSeriesContainer getDailyEvapTimeSeries() {
        TimeSeriesContainer evapTs = this.getHourlyEvapTimeSeries();
        TimeSeriesMath tsAcc = null;

        try {
            TimeSeriesMath tsMath = new TimeSeriesMath(evapTs);
            tsAcc = (TimeSeriesMath)tsMath.transformTimeSeries("1DAY", "", "ACC", false);
        } catch (HecMathException var4) {
            HecMathException hme = var4;
            LOGGER.log(Level.SEVERE, "HecMathException occurred while converting hourly total evap to daily total evap", hme);
        }

        return tsAcc != null ? tsAcc.getContainer() : null;
    }

    public TimeSeriesContainer getDailyEvapFlowTimeSeries() {
        TimeSeriesContainer dailyEvapTs = this.getDailyEvapTimeSeries();
        TimeSeriesContainer dailyEvapFlowTs = (TimeSeriesContainer)dailyEvapTs.clone();
        double evap_to_meters = 1.0;

        try {
            evap_to_meters = Units.convertUnits(1.0, dailyEvapTs.units.toString(), "m");
        } catch (UnitsConversionException var22) {
            UnitsConversionException ue = var22;
            LOGGER.log(Level.SEVERE, "Unable to convert " + dailyEvapTs.units + " to m", ue);
            return null;
        }

        double undef = -3.4028234663852886E38;
        double lastValidElev = undef;
        HecTime hecTime = new HecTime(1);

        for(int i = 0; i < dailyEvapFlowTs.numberValues; ++i) {
            int itime = dailyEvapFlowTs.times[i];
            hecTime.set(itime);
            double dailyEvap = dailyEvapTs.values[i];
            if (RMAConst.isValidValue(dailyEvap)) {
                double elev = this._reservoir.getCurrentElevation(hecTime);
                double currentElev;
                if (RMAConst.isValidValue(elev)) {
                    lastValidElev = elev;
                    currentElev = elev;
                } else {
                    if (!RMAConst.isValidValue(lastValidElev)) {
                        dailyEvapFlowTs.values[i] = undef;
                        continue;
                    }

                    currentElev = lastValidElev;
                }

                double areaMetersSq = this._reservoir.intArea(currentElev);
                double dailyEvapFlow = areaMetersSq * dailyEvap * evap_to_meters;
                dailyEvapFlowTs.values[i] = dailyEvapFlow / 86400.0;
            } else {
                dailyEvapFlowTs.values[i] = undef;
            }
        }

        dailyEvapFlowTs.units = "cms";
        dailyEvapFlowTs.parameter = "FLOW-EVAP";
        DSSPathString dsspath = new DSSPathString(dailyEvapFlowTs.fullName);
        dsspath.setCPart(dailyEvapFlowTs.parameter);
        dailyEvapFlowTs.fullName = dsspath.getPathname();
        dailyEvapFlowTs.type = "Ave";
        return dailyEvapFlowTs;
    }

    public TimeSeriesContainer[] getTemperatureProfileTs(double surfaceDepth, double bottomDepth, double intervalDepth) {
        if (intervalDepth < 0.01) {
            return null;
        } else if (this._wtempProfiles != null && this._wtempProfiles.length >= 1) {
            if (bottomDepth < 0.0) {
                bottomDepth = -bottomDepth;
            }

            double rprofs = (bottomDepth - surfaceDepth) / intervalDepth;
            int nprofs = (int)(rprofs + 0.01) + 1;
            TimeSeriesContainer wtTemplate = (TimeSeriesContainer)this._surfaceTempTsc.clone();
            DSSPathString dsspath = new DSSPathString(this._surfaceTempTsc.fullName);
            dsspath.setCPart("Temp-Water");
            wtTemplate.parameter = "Temp-Water";
            wtTemplate.fullName = dsspath.getPathname();
            TimeSeriesContainer[] wtprofTs = new TimeSeriesContainer[nprofs];
            double dep = surfaceDepth;

            for(int i = 0; i < nprofs; ++i) {
                wtprofTs[i] = (TimeSeriesContainer)wtTemplate.clone();
                int iint = (int)dep;
                double frac = dep % 1.0;
                int ifrac = (int)Math.round(frac * 10.0);
                String fstr = Integer.toString(iint) + "," + Integer.toString(ifrac);
                String paramstr = wtTemplate.parameter + "-" + fstr + "m";
                wtprofTs[i].parameter = paramstr;
                dsspath.setCPart(paramstr);
                wtprofTs[i].fullName = dsspath.getPathname();
                dep += intervalDepth;
            }

            double undef = -3.4028234663852886E38;

            int nvals;
            for(nvals = 0; nvals < nprofs; ++nvals) {
                Arrays.fill(wtprofTs[nvals].values, undef);
            }

            nvals = this._wtempProfiles.length;
            int layers = this._wtempProfiles[0].length;
            int resj = this._reservoir._resj;

            for(int i = resj; i >= 0; --i) {
                int ilayer = resj - i;

                for(int itime = 0; itime < nvals; ++itime) {
                    double tval = this._wtempProfiles[itime][i];
                    wtprofTs[ilayer].values[itime] = tval;
                }
            }

            return wtprofTs;
        } else {
            return null;
        }
    }

    public TimeSeriesContainer[] getDailyTemperatureProfileTs(double surfaceDepth, double bottomDepth, double intervalDepth) {
        LOGGER.log(Level.SEVERE, "getDailyTemperatureProfileTs");
        TimeSeriesContainer[] hourlyTsArray = this.getTemperatureProfileTs(surfaceDepth, bottomDepth, intervalDepth);
        int nlayers = hourlyTsArray.length;
        TimeSeriesContainer[] dayTsArray = new TimeSeriesContainer[nlayers];
        int icnt = 0;
        TimeSeriesContainer[] var11 = hourlyTsArray;
        int var12 = hourlyTsArray.length;

        for(int var13 = 0; var13 < var12; ++var13) {
            TimeSeriesContainer hourlyTsc = var11[var13];
            if (Arrays.stream(hourlyTsc.values).noneMatch(RMAConst::isValidValue)) {
                LOGGER.log(Level.FINE, () -> {
                    return "No data found for " + hourlyTsc.parameter;
                });
                ++icnt;
            } else {
                try {
                    TimeSeriesMath tsMath = new TimeSeriesMath(hourlyTsc);
                    tsMath.getContainer().type = "INST-VAL";
                    TimeSeriesMath tsMath24 = (TimeSeriesMath)tsMath.transformTimeSeries("1DAY", "", "INT", false);
                    TimeSeriesContainer dailyTsc = tsMath24.getContainer();
                    this.verifyHourlyAndDailyMatch(hourlyTsc, dailyTsc);
                    dayTsArray[icnt] = dailyTsc;
                } catch (HecMathException var18) {
                    HecMathException ex = var18;
                    LOGGER.log(Level.SEVERE, ex, () -> {
                        return "Exception occurred while transforming hourly temperature profile data to daily for " + hourlyTsc.parameter + ".";
                    });
                    dayTsArray[icnt] = null;
                }

                ++icnt;
            }
        }

        return dayTsArray;
    }

    private void verifyHourlyAndDailyMatch(TimeSeriesContainer hourlyTsc, TimeSeriesContainer dailyTsc) {
        if (LOGGER.isLoggable(Level.FINE)) {
            LOGGER.log(Level.SEVERE, () -> {
                return "Checking Daily and Hourly timeseries " + dailyTsc.fullName + " for discrepancies in values and times";
            });
            List<Integer> hourlyTimes = (List)Arrays.stream(hourlyTsc.times).boxed().collect(Collectors.toList());
            int dailyIndex = 0;
            int[] var5 = dailyTsc.times;
            int var6 = var5.length;

            for(int var7 = 0; var7 < var6; ++var7) {
                int dailyTime = var5[var7];
                HecTime dailyHecTime = dailyTsc.getTimes().elementAt(dailyIndex);
                int hourlyIndex = hourlyTimes.indexOf(dailyTime);
                if (hourlyIndex == -1) {
                    LOGGER.log(Level.SEVERE, () -> {
                        return "Hourly data doesn't contain a time for " + dailyHecTime.dateAndTime();
                    });
                } else {
                    double hourlyValue = hourlyTsc.values[hourlyIndex];
                    double dailyValue = dailyTsc.values[dailyIndex];
                    if (!RmaMath.equals(hourlyValue, dailyValue, 9.999999747378752E-6)) {
                        LOGGER.log(Level.SEVERE, () -> {
                            return "Hourly and Daily value don't match at " + dailyHecTime.dateAndTime() + System.lineSeparator() + "\tExpected: " + hourlyValue + System.lineSeparator() + "\tReceived: " + dailyValue;
                        });
                    }

                    ++dailyIndex;
                }
            }
        }

    }

    public static TimeSeriesContainer[] generateDailyTemperatureProfileTs(String locationName, String versionName, String startDate, String endDate, double tval, double surfaceDepth, double bottomDepth, double intervalDepth) {
        if (intervalDepth < 0.01) {
            return null;
        } else {
            if (bottomDepth < 0.0) {
                bottomDepth = -bottomDepth;
            }

            double rprofs = (bottomDepth - surfaceDepth) / intervalDepth;
            int nprofs = (int)(rprofs + 0.01) + 1;
            int intervalMinutes = 1440;
            int offsetMinutes = 0;
            TimeSeriesContainer baseTsc = null;

            try {
                TimeSeriesMath tsMath = (TimeSeriesMath)TimeSeriesMath.generateRegularIntervalTimeSeries(startDate, endDate, intervalMinutes, offsetMinutes, tval);
                baseTsc = tsMath.getContainer();
            } catch (HecMathException var30) {
                HecMathException hme = var30;
                LOGGER.log(Level.SEVERE, "Exception occurred while generating daily temperature profile data.", hme);
                return null;
            }

            baseTsc.location = locationName;
            baseTsc.version = versionName;
            baseTsc.interval = intervalMinutes;
            baseTsc.units = "C";
            baseTsc.type = "Inst";
            baseTsc.units = "C";
            DSSPathString dsspath = new DSSPathString("", locationName, "", "", "1DAY", versionName);
            dsspath.setCPart("Temp-Water");
            baseTsc.parameter = "Temp-Water";
            baseTsc.fullName = dsspath.getPathname();
            TimeSeriesContainer wtTemplate = (TimeSeriesContainer)baseTsc.clone();
            TimeSeriesContainer[] wtprofTs = new TimeSeriesContainer[nprofs];
            double dep = surfaceDepth;

            for(int i = 0; i < nprofs; ++i) {
                wtprofTs[i] = (TimeSeriesContainer)wtTemplate.clone();
                int iint = (int)dep;
                double frac = dep % 1.0;
                int ifrac = (int)Math.round(frac * 10.0);
                String fstr = Integer.toString(iint) + "," + Integer.toString(ifrac);
                String paramstr = wtTemplate.parameter + "-" + fstr + "m";
                wtprofTs[i].parameter = paramstr;
                dsspath.setCPart(paramstr);
                wtprofTs[i].fullName = dsspath.getPathname();
                dep += intervalDepth;
            }

            return wtprofTs;
        }
    }
}
