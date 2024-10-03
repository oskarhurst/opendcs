package decodes.cwms.newresevap;

public enum WindShearMethod
{
    DONELAN("Donelan"),
    FISCHER("Fischer");

    private final String _name;

    private WindShearMethod(String name) {
        this._name = name;
    }

    public String toString() {
        return this._name;
    }

    public static WindShearMethod fromString(String name) {
        WindShearMethod[] var1 = values();
        int var2 = var1.length;

        for(int var3 = 0; var3 < var2; ++var3) {
            WindShearMethod windShearMethod = var1[var3];
            if (windShearMethod._name.equalsIgnoreCase(name)) {
                return windShearMethod;
            }
        }

        throw new IllegalArgumentException("No wind shear method with name " + name + " found");
    }
}
