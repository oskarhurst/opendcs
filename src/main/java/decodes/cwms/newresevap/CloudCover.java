package decodes.cwms.newresevap;

public class CloudCover {
    public double fractionCloudCover;
    public double height;
    public CloudHeightType cloudType;
    private static final double MISSING_VALUE_LOW = 0.54;
    private static final double MISSING_VALUE_MED = 0.0;
    private static final double MISSING_VALUE_HI = 0.0;

    public CloudCover(double fraction, double height, CloudHeightType type) {
        this.fractionCloudCover = fraction;
        this.height = height;
        this.cloudType = type;
    }

    public int getCloudTypeFlag() {
        int[] var10000 = new int[]{2, 3, 4};
        if (this.cloudType == CloudCover.CloudHeightType.height_low) {
            return 4;
        } else if (this.cloudType == CloudCover.CloudHeightType.height_med) {
            return 3;
        } else {
            return this.cloudType == CloudCover.CloudHeightType.height_high ? 2 : 3;
        }
    }

    public String getTypeName() {
        if (this.cloudType == CloudCover.CloudHeightType.height_low) {
            return "Height Low";
        } else if (this.cloudType == CloudCover.CloudHeightType.height_med) {
            return "Height Med";
        } else {
            return this.cloudType == CloudCover.CloudHeightType.height_high ? "Height High" : "";
        }
    }

    public double getDefaultFractionCloudCover() {
        if (this.cloudType == CloudCover.CloudHeightType.height_low) {
            return 0.54;
        } else if (this.cloudType == CloudCover.CloudHeightType.height_med) {
            return 0.0;
        } else {
            return this.cloudType == CloudCover.CloudHeightType.height_high ? 0.0 : 0.0;
        }
    }

    public static enum CloudHeightType {
        height_low,
        height_med,
        height_high;

        private CloudHeightType() {
        }
    }
}