/**
 *
 */
package org.theseed.reports;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.theseed.genome.Feature;
import org.theseed.genome.Genome;
import org.theseed.io.LineReader;

/**
 * @author Bruce Parrello
 *
 */
public abstract class ApplyKmerReporter extends BaseReporter {

    // FIELDS
    /** map of role IDs to output column indices */
    private Map<String, Integer> roleIdxMap;

    /**
     * Construct a report object.
     *
     * @param output	output stream for the report
     */
    public ApplyKmerReporter(OutputStream output) {
        super(output);
    }

    /**
     * Initialize the report.
     *
     * @param rolesToUse	file containing the interesting role IDs, in order, in the first column
     *
     * @throws IOException
     */
    public void initReport(File rolesToUse) throws IOException {
        this.roleIdxMap = new HashMap<String, Integer>(2500);
        try (LineReader roleStream = new LineReader(rolesToUse)) {
            int idx = 1;
            for (String line : roleStream) {
                String role = StringUtils.substringBefore(line, "\t");
                this.roleIdxMap.put(role, idx);
                idx++;
            }
        }
        this.openReport();
    }

    /**
     * Allow the subclass to initialize the report.
     */
    protected abstract void openReport();

    /**
     * Begin reporting on a genome.
     *
     * @param genome	genome of interest
     */
    public abstract void openGenome(Genome genome);

    /**
     * Record hits on a feature.
     *
     * @param feat		feature hit
     * @param role		ID of the role hit
     * @param count		number of hits
     */
    public abstract void recordFeature(Feature feat, String role, int count);

    /**
     * End processing on a genome.
     */
    public abstract void closeGenome();

    /**
     * Finish the report.
     */
    public abstract void closeReport();

    /**
     * @return the column index for the specified role ID, or 0 if the role is not interesting
     *
     * @param roleId	ID of the indicated role
     */
    public int getRoleIdx(String roleId) {
        Integer retVal = this.roleIdxMap.getOrDefault(roleId, 0);
        return (int) retVal;
    }

    /**
     * @return the number of roles in the map
     */
    public int getNumRoles() {
        return this.roleIdxMap.size();
    }

    /**
     * Enum for types of reports.
     */
    public static enum Type {
        VERIFY, APPLY;

        /**
         * Create a reporting object of the specified type.
         */
        public ApplyKmerReporter create(OutputStream output) {
            ApplyKmerReporter retVal = null;
            switch (this) {
            case APPLY:
                retVal = new DefaultApplyKmerReporter(output);
                break;
            case VERIFY:
                retVal = new VerifyApplyKmerReporter(output);
                break;
            }
            return retVal;
        }
    }
}
