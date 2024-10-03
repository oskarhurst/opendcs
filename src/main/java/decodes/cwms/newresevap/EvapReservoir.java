/* 
 * Copyright (c) 2018
 * United States Army Corps of Engineers - Hydrologic Engineering Center (USACE/HEC)
 * All Rights Reserved.  USACE PROPRIETARY/CONFIDENTIAL.
 * Source may not be released without written approval from HEC
 */
package decodes.cwms.newresevap;

import decodes.cwms.HecConstants;
import decodes.cwms.resevapcalc.ReservoirLocationInfo;
import decodes.cwms.resevapcalc.WindShearMethod;
import decodes.tsdb.CTimeSeries;
import decodes.util.DecodesException;
import hec.data.RatingException;
import hec.data.cwmsRating.RatingSet;

import java.io.*;
import java.sql.Connection;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

//import rma.util.RMAConst;

/**
 * This class holds reservoir specific values (lat, lon, elev-area)
 * and the layers info for computing the
 * water temperature profile.
 * 
 */
public class EvapReservoir
{
    private static final Logger LOGGER = Logger.getLogger(EvapReservoir.class.getName());
    public static final int MAX_RES_AREAS = 50000;
    public static final int NLAYERS = 1000;
    public static final double FT_TO_M = 0.3048;
    public static final double ACRES_TO_M2 = 4046.8564224;
    public int _numberLayers;
    public double _secchiDepth;
    public double _attenuationConst;
    public double _zeroElevaton;
    public double _lat;
    public double _long;
    public WindShearMethod _windShearMethod;
    public double _thermalDiffusivityCoefficient;
    protected double _ru;
    protected double _rt;
    protected double _rq;
    double _gmtOffset;
    String _tzString;
    public double _elev;
    public double _instrumentHeight;
    public double[] _elevA;
    public double[] _surfA;
    public int _nSurfArea;
    public int _nResBottom;
    public double _depth;
    public double _surfArea;
    public double[] _zd;
    public double[] _delz;
    public double[] _ztop;
    public double[] _zarea;
    public double[] _zelev;
    public double[] _zvol;
    public double[] _rhow;
    public double[] _cp;
    public double[] _wt;
    public double[] _kz;
    public TimeSeriesContainer _elevationTsc;
    int _resj;
    int _resj_old;
    int _resj1;
    boolean _isEnglish;
    boolean _isDewpoint;
    boolean _elevTS_flag;
    String _name;
    BufferedWriter _debugout;

    public EvapReservoir() {
        this._windShearMethod = WindShearMethod.DONELAN;
        this._thermalDiffusivityCoefficient = 1.2;
        this._zd = new double[1000];
        this._delz = new double[1000];
        this._ztop = new double[1000];
        this._zarea = new double[1000];
        this._zelev = new double[1000];
        this._zvol = new double[1000];
        this._rhow = new double[1000];
        this._cp = new double[1000];
        this._wt = new double[1000];
        this._kz = new double[1000];
        this._isEnglish = true;
        this._isDewpoint = false;
        this._elevTS_flag = false;
        this._debugout = null;
    }

    public void setInputDataIsEnglish(boolean tf) {
        this._isEnglish = tf;
    }

    public void setName(String name) {
        this._name = name;
    }

    public String getName() {
        return this._name;
    }

    public void setTimeZoneString(String tzStr) {
        this._tzString = tzStr;
    }

    public void setLatLon(double lat, double lon) {
        this._lat = lat;
        this._long = lon;
    }

    public void setGmtOffset(double gmtOffset) {
        this._gmtOffset = gmtOffset;
    }

    public double getGmtOffset() {
        return this._gmtOffset;
    }

    public void setDebugFile(BufferedWriter debugout) {
        this._debugout = debugout;
    }

    public double getCurrentElevation(HecTime hecTime) {
        if (this._elevationTsc == null) {
            return this._elev;
        } else {
            double elev = this._elevationTsc.getValue(hecTime);
            if (!RMAConst.isValidValue(elev)) {
                return elev;
            } else {
                try {
                    elev = Units.convertUnits(elev, this._elevationTsc.units.toString(), "m");
                } catch (UnitsConversionException var5) {
                    UnitsConversionException ue = var5;
                    LOGGER.log(Level.SEVERE, "Exception occurred while transforming hourly temperature profile data to daily for .", ue);
                }

                return elev;
            }
        }
    }

    public void setInstrumentHeights(double windHeight, double tempHeigth, double relHeight) {
        double sclfct = 1.0;
        if (this._isEnglish) {
            sclfct = 0.3048;
        }

        this._ru = windHeight * sclfct;
        this._rt = tempHeigth * sclfct;
        this._rq = relHeight * sclfct;
        this._instrumentHeight = windHeight * sclfct;
    }

    public void setSecchi(double secchi) {
        this._secchiDepth = secchi;
        if (this._isEnglish) {
            this._secchiDepth *= 0.3048;
        }

        this._attenuationConst = 1.7 / this._secchiDepth;
    }

    public void setWindShearMethod(WindShearMethod method) {
        this._windShearMethod = method;
    }

    public void setThermalDiffusivityCoefficient(double thermalDiffusivityCoefficient) {
        this._thermalDiffusivityCoefficient = thermalDiffusivityCoefficient;
    }

    public boolean readTextFile(File textFile) {
        List<String> inputLines = new ArrayList();
        if (textFile != null && textFile.exists()) {
            int icnt = 0;

            try {
                BufferedReader in = new BufferedReader(new FileReader(textFile));
                Throwable var5 = null;

                try {
                    for(String line = in.readLine(); line != null; line = in.readLine()) {
                        ++icnt;
                        inputLines.add(line);
                    }
                } catch (Throwable var18) {
                    var5 = var18;
                    throw var18;
                } finally {
                    if (in != null) {
                        if (var5 != null) {
                            try {
                                in.close();
                            } catch (Throwable var17) {
                                var5.addSuppressed(var17);
                            }
                        } else {
                            in.close();
                        }
                    }

                }
            } catch (IOException var20) {
                IOException e = var20;
                LOGGER.log(Level.SEVERE, "Exception:  Error Reading bc file ", e);
                return false;
            }

            int i;
            for(i = 0; i < inputLines.size(); ++i) {
                String line = (String)inputLines.get(i);
                if (line.trim().length() > 1) {
                    String[] parts = line.split(":");
                    if (parts != null && parts.length == 2) {
                        if (parts[0].toUpperCase().contains("RESERVOIR")) {
                            this._name = parts[1];
                        } else if (parts[0].toUpperCase().contains("METRIC")) {
                            this._isEnglish = false;
                        } else if (parts[0].toUpperCase().contains("DEWPOINT")) {
                            this._isDewpoint = true;
                        } else if (parts[0].toUpperCase().contains("SECCHI")) {
                            this._secchiDepth = Double.parseDouble(parts[1]);
                            if (this._isEnglish) {
                                this._secchiDepth *= 0.3048;
                            }

                            this._attenuationConst = 1.7 / this._secchiDepth;
                        } else if (parts[0].toUpperCase().contains("ZERO ELEV")) {
                            this._zeroElevaton = Double.parseDouble(parts[1]);
                        } else if (parts[0].toUpperCase().contains("LAT")) {
                            this._lat = Double.parseDouble(parts[1]);
                        } else {
                            double npairs;
                            if (parts[0].toUpperCase().contains("LONG")) {
                                npairs = Double.parseDouble(parts[1]);
                                this._long = Math.abs(npairs);
                            } else if (parts[0].toUpperCase().contains("GMT")) {
                                this._gmtOffset = Double.parseDouble(parts[1]);
                            } else if (parts[0].toUpperCase().contains("ELEV")) {
                                this._elev = Double.parseDouble(parts[1]);
                            } else if (parts[0].toUpperCase().contains("INSTRUMENT HEIGHT")) {
                                this._instrumentHeight = Double.parseDouble(parts[1]);
                                this._ru = this._instrumentHeight;
                                this._rt = this._instrumentHeight;
                                this._rq = this._instrumentHeight;
                            } else if (parts[0].toUpperCase().contains("SURFACE AREA")) {
                                npairs = Double.parseDouble(parts[1]);
                                this._nSurfArea = (int)npairs;
                                if (this._nSurfArea > 50000) {
                                    return false;
                                }

                                double[] resElevations = new double[this._nSurfArea];
                                double[] resAreas = new double[this._nSurfArea];
                                ++i;
                                line = (String)inputLines.get(i);
                                parts = line.split(",");
                                if (parts.length != this._nSurfArea) {
                                    return false;
                                }

                                int j;
                                for(j = 0; j < parts.length; ++j) {
                                    resElevations[j] = Double.parseDouble(parts[j]);
                                }

                                ++i;
                                line = (String)inputLines.get(i);
                                parts = line.split(",");
                                if (parts.length != this._nSurfArea) {
                                    return false;
                                }

                                for(j = 0; j < parts.length; ++j) {
                                    resAreas[j] = Double.parseDouble(parts[j]);
                                }

                                this._elevA = resElevations;
                                this._surfA = resAreas;
                                break;
                            }
                        }
                    }
                }
            }

            if (this._isEnglish) {
                i = this._elevA.length;

                for(int i = 0; i < i; ++i) {
                    double[] var10000 = this._elevA;
                    var10000[i] *= 0.3048;
                    var10000 = this._surfA;
                    var10000[i] *= 4046.8564224;
                }

                this._ru *= 0.3048;
                this._rt *= 0.3048;
                this._rq *= 0.3048;
                this._instrumentHeight *= 0.3048;
                if (this._zeroElevaton > -1.0) {
                    this._zeroElevaton *= 0.3048;
                }

                this._elev *= 0.3048;
            }

            double surfArea = this.intArea(this._elev);
            this._surfArea = surfArea / 1000000.0;
            return true;
        } else {
            return false;
        }
    }

    public ReservoirLocationInfo getReservoirLocationInfo() {
        ReservoirLocationInfo resInfo = new ReservoirLocationInfo();
        resInfo.lat = this._lat;
        resInfo.lon = this._long;
        resInfo.instrumentHeight = this._instrumentHeight;
        resInfo.ru = this._ru;
        resInfo.rt = this._rt;
        resInfo.rq = this._rq;
        resInfo.gmtOffset = this._gmtOffset;
        return resInfo;
    }

    public boolean readInitTemperature(File textFile) {
        List<String> inputLines = new ArrayList();
        if (textFile != null && textFile.exists()) {
            int icnt = 0;

            String line;
            try {
                BufferedReader in = new BufferedReader(new FileReader(textFile));
                Throwable var5 = null;

                try {
                    for(line = in.readLine(); line != null; line = in.readLine()) {
                        ++icnt;
                        inputLines.add(line);
                    }
                } catch (Throwable var15) {
                    var5 = var15;
                    throw var15;
                } finally {
                    if (in != null) {
                        if (var5 != null) {
                            try {
                                in.close();
                            } catch (Throwable var14) {
                                var5.addSuppressed(var14);
                            }
                        } else {
                            in.close();
                        }
                    }

                }
            } catch (IOException var17) {
                IOException e = var17;
                LOGGER.log(Level.SEVERE, "Exception:  Error Reading intial water temperature file", e);
                return false;
            }

            icnt = 0;
            double[] wt = new double[1000];

            for(int i = 0; i < inputLines.size(); ++i) {
                line = (String)inputLines.get(i);
                if (line.trim().length() > 1) {
                    wt[icnt] = Double.parseDouble(line);
                    ++icnt;
                }
            }

            this._wt = wt;
            return true;
        } else {
            return false;
        }
    }

    public boolean setInitWaterTemperature(double waterTemp, int resj) {
        for(int i = 0; i <= resj; ++i) {
            this._wt[i] = waterTemp;
        }

        return true;
    }

    public boolean setInitWaterTemperatureProfile(double[] wt, int resj) {
        for(int i = 0; i <= resj; ++i) {
            this._wt[i] = wt[i];
        }

        return true;
    }

    public boolean setElevAreaCurve(double[] elev, double[] area, int nvals, boolean isEnglish) {
        this._elevA = new double['썐'];
        this._surfA = new double['썐'];
        Arrays.fill(this._elevA, 0.0);
        Arrays.fill(this._surfA, 0.0);
        double elevSclfct = 1.0;
        double areaSclfct = 1.0;
        if (isEnglish) {
            elevSclfct = 0.3048;
            areaSclfct = 4046.8564224;
        }

        for(int i = 0; i < nvals; ++i) {
            this._elevA[i] = elev[i] * elevSclfct;
            this._surfA[i] = area[i] * areaSclfct;
        }

        this._nSurfArea = nvals;
        return true;
    }

    public boolean setElevation(double elev) {
        this._elev = elev;
        if (this._isEnglish) {
            this._elev *= 0.3048;
        }

        double surfArea = this.intArea(this._elev);
        this._surfArea = surfArea / 1000000.0;
        return true;
    }

    public boolean setElevationMeters(double elev) {
        this._elev = elev;
        return true;
    }

    public double getElevation() {
        return this._elev;
    }

    public boolean setZeroElevation(double elev) {
        this._zeroElevaton = elev;
        if (this._isEnglish && elev > -1.0) {
            this._zeroElevaton *= 0.3048;
        }

        return true;
    }

    public void setElevationTs(TimeSeriesContainer tsc) {
        this._elevationTsc = tsc;
    }

    public boolean initRes() {
        double wsel = this._elev;
        int nResBottom = -1;

        for(int i = 0; i < this._nSurfArea; ++i) {
            LOGGER.log(Level.FINE, "i, surfA, nResBottom  " + i + ", " + nResBottom + ", " + this._surfA[i]);
            if (this._surfA[i] > 0.0) {
                if (i == 0) {
                    nResBottom = 0;
                } else {
                    nResBottom = i - 1;
                }

                if (nResBottom > -1) {
                    break;
                }
            }
        }

        double depth;
        if (this._zeroElevaton == -1.0) {
            depth = wsel - this._elevA[nResBottom];
        } else {
            int nResTest = this.locate(this._zeroElevaton);
            if (nResTest > nResBottom) {
                nResBottom = nResTest;
            } else {
                this._zeroElevaton = this._elevA[nResBottom];
            }

            depth = wsel - this._zeroElevaton;
        }

        if (depth <= 0.0) {
            String msg = "Water Elev less than Reservoir Bottom Elevation";
            LOGGER.log(Level.WARNING, msg);
            return false;
        } else {
            this._depth = depth;
            this._nResBottom = nResBottom;
            LOGGER.log(Level.FINE, "nResBottom  {0}, {1}", new Object[]{this._nResBottom, this._nSurfArea});
            double surfArea = this.intArea(this._elev);
            this._surfArea = surfArea / 1000000.0;
            return true;
        }
    }

    public int getResj() {
        return this._resj;
    }

    public boolean resSetup() {
        return this.resSetup(false);
    }

    public boolean resSetup(boolean updateDepth) {
        double[] zdx = new double[1000];
        double[] delzx = new double[1000];
        double[] ztop_depth = new double[1000];
        double[] zelevx = new double[1000];
        double[] zareax = new double[1000];
        double[] zvolx = new double[1000];
        double wsel = this._elev;
        if (updateDepth) {
            this._depth = wsel - this._elevA[this._nResBottom];
        }

        double zdepth = 0.25;

        int j;
        for(j = 0; zdepth <= this._depth; ++j) {
            zdx[j] = zdepth;
            delzx[j] = 0.5;
            zdepth += 0.5;
            ztop_depth[j] = (double)j * 0.5;
            zelevx[j] = wsel - (double)j * 0.5;
            zareax[j] = this.intArea(zelevx[j]);
        }

        int resj = j - 1;
        LOGGER.log(Level.FINE, " resj,  depth, wsel {0}  {1}  {2}", new Object[]{resj, this._depth, wsel});
        if (zdx[resj] < this._depth) {
            zdx[resj] = this._depth;
            delzx[resj] = zdx[resj] - zdx[resj - 1] - 0.5 * delzx[resj - 1];
            ztop_depth[resj] = this._depth - delzx[resj];
            zelevx[resj] = this._elevA[this._nResBottom] + delzx[resj];
            zareax[resj] = this.intArea(zelevx[resj]);
        }

        for(j = 0; j < resj; ++j) {
            zvolx[j] = 0.5 * (zareax[j] + zareax[j + 1]) * delzx[j];
        }

        zvolx[resj] = 0.5 * zareax[resj] * delzx[resj];
        int i = 0;

        for(j = resj; j >= 0; --j) {
            this._zd[i] = zdx[j];
            this._delz[i] = delzx[j];
            this._ztop[i] = this._depth - ztop_depth[j];
            this._zarea[i] = zareax[j];
            this._zelev[i] = zelevx[j];
            this._zvol[i] = zvolx[j];
            ++i;
        }

        this._resj_old = this._resj;
        this._resj = resj;
        if (this._debugout != null) {
            this.writeLayerDebugInfo();
        }

        return true;
    }

    public void writeLayerDebugInfo() {
        if (this._debugout != null) {
            NumberFormat nf10_3 = new NumberFormat("%10.3f");
            NumberFormat nf13_4 = new NumberFormat("%13.4f");
            NumberFormat nfi_3 = new NumberFormat("%3d");

            try {
                String depthstr = "resj+1 : " + this._resj + 1 + "      wsel_ft" + nf10_3.form(this._elev / 0.3048) + "      wsel" + nf10_3.form(this._elev) + "     depth" + nf10_3.form(this._depth);
                this._debugout.write(depthstr);
                this._debugout.newLine();
                String heading = " j               zd            delz            ztop           zarea           zelev            zvol";
                StringBuffer strbuf = new StringBuffer();
                this._debugout.write(heading);
                this._debugout.newLine();

                for(int j = 0; j <= this._resj; ++j) {
                    strbuf.setLength(0);
                    strbuf.append(nfi_3.form((long)(j + 1)));
                    strbuf.append(nf13_4.form(this._zd[j]) + "  ");
                    strbuf.append(nf13_4.form(this._delz[j]) + "  ");
                    strbuf.append(nf13_4.form(this._ztop[j]) + "  ");
                    strbuf.append(nf13_4.form(this._zarea[j]) + "  ");
                    strbuf.append(nf13_4.form(this._zelev[j]) + "  ");
                    strbuf.append(nf13_4.form(this._zvol[j]) + "  ");
                    this._debugout.write(strbuf.toString());
                    this._debugout.newLine();
                }
            } catch (IOException var8) {
                IOException ioe = var8;
                LOGGER.log(Level.FINE, "Exception occurred while writing layer debug info.", ioe);
            }

        }
    }

    public int locate(double x) {
        int n = this._nSurfArea - 1;
        int jl = -1;
        int ju = n + 1;

        while(true) {
            while(ju - jl > 1) {
                int jm = (ju + jl) / 2;
                if (this._elevA[n] >= this._elevA[0] && x >= this._elevA[jm]) {
                    jl = jm;
                } else {
                    ju = jm;
                }
            }

            int j;
            if (x == this._elevA[0]) {
                j = 0;
            } else if (x == this._elevA[n]) {
                j = n - 1;
            } else {
                j = jl;
            }

            return j;
        }
    }

    public double intArea(double el) {
        int j = this.locate(el);
        double w = this._surfA[j] + (el - this._elevA[j]) / (this._elevA[j + 1] - this._elevA[j]) * (this._surfA[j + 1] - this._surfA[j]);
        return w;
    }

}
