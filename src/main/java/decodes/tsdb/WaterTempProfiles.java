package decodes.tsdb;

import decodes.db.Constants;
import decodes.db.Site;
import decodes.db.SiteName;
import ilex.var.TimedVariable;
import opendcs.dai.TimeSeriesDAI;

import java.util.ArrayList;
import java.util.Date;

public class WaterTempProfiles extends DataCollection {

    private TimeSeriesDAI timeSeriesDAO = null;

    /** The time series */
    private ArrayList<CTimeSeries> tseries;

    /** Handle storing tasklist Record Ranges for computation processor. */
    private RecordRangeHandle rrHandle;

    private double startDepth;
    private double increment;

    /** Constructor -- builds an empty collection with a null handle. */
    public WaterTempProfiles(TimeSeriesDAI DAO, double start, double incr){
        tseries = new ArrayList<CTimeSeries>();
        rrHandle = null;
        timeSeriesDAO = DAO;
        startDepth = start;
        increment = incr;
    }
    public WaterTempProfiles(ArrayList<CTimeSeries> profiles, TimeSeriesDAI DAO, double start, double incr){
        tseries = profiles;
        rrHandle = null;
        timeSeriesDAO = DAO;
        startDepth = start;
        increment = incr;
    }

    public WaterTempProfiles(TimeSeriesDAI DAO, String wtpId, Date since, Date until, double start, double incr) throws DbIoException, NoSuchObjectException {
        tseries = new ArrayList<CTimeSeries>();
        rrHandle = null;
        timeSeriesDAO = DAO;
        startDepth = start;
        increment = incr;
        boolean loading = true;
        double currentDepth = startDepth;
        TimeSeriesIdentifier tsid;
        tsid = timeSeriesDAO.getTimeSeriesIdentifier(wtpId);
        while(loading){
            try {
                TimeSeriesIdentifier newtsid = tsid.copyNoKey();
                Site newsite =  new Site();
                newsite.copyFrom(tsid.getSite());
                SiteName strsite = newsite.getName(Constants.snt_CWMS);
                strsite.setNameValue(strsite.getNameValue()+currentDepth);
                newtsid.setSite(newsite);
                CTimeSeries cts = timeSeriesDAO.makeTimeSeries(newtsid);
                int n = timeSeriesDAO.fillTimeSeries(cts, since, until);
                if (n == 0){
                    loading = false;
                }
                else{
                    tseries.add(cts);
                }
                currentDepth += increment;
            }
            catch (BadTimeSeriesException ex)
            {
                loading = false;
            }
        }

    }

    public void SaveProfiles(){
        for (CTimeSeries tsery : tseries) {
            try {
                timeSeriesDAO.saveTimeSeries(tsery);
            } catch (Exception ex) {
                String msg = "Error saving water temperature profile data: " + ex;
                //warning(msg);
                System.err.print(msg);
                ex.printStackTrace(System.err);
            }
        }
    }



}