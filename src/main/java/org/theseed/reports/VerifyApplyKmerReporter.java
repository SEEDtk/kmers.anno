/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;

import org.theseed.genome.Feature;
import org.theseed.genome.Genome;

/**
 * This report compares the roles found to the functions already present in the genome.
 *
 * @author Bruce Parrello
 *
 */
public class VerifyApplyKmerReporter extends ApplyKmerReporter {

    // FIELDS
    /** ID of the current genome */
    private String genomeId;

    /**
     * Attach the specified output stream to this report.
     *
     * @param output	output stream for the report
     */
    public VerifyApplyKmerReporter(OutputStream output) {
        super(output);
    }

    @Override
    public void openReport() {
        this.println("genome_id\tpeg_id\trole\thits\tfunction");
    }

    @Override
    public void openGenome(Genome genome) {
        this.genomeId = genome.getId();
    }

    @Override
    public void recordFeature(Feature feat, String role, int count) {
        this.print("%s\t%s\t%s\t%d\t%s", this.genomeId, feat.getId(), role, count, feat.getFunction());
    }

    @Override
    public void closeGenome() {
    }

    @Override
    public void closeReport() {
    }

}
