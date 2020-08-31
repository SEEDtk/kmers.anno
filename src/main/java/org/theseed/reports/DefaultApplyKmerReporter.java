/**
 *
 */
package org.theseed.reports;

import java.io.OutputStream;
import java.util.Arrays;
import java.util.stream.Collectors;

import org.theseed.genome.Feature;
import org.theseed.genome.Genome;

/**
 * @author Bruce Parrello
 *
 */
public class DefaultApplyKmerReporter extends ApplyKmerReporter {

    // FIELDS
    /** role counts for the current genome */
    private int[] roleCounts;
    /** ID of the current genome */
    private String genomeId;

    /**
     * @param output
     */
    public DefaultApplyKmerReporter(OutputStream output) {
        super(output);
    }

    @Override
    public void openReport() {
        this.roleCounts = new int[this.getNumRoles()];
    }

    @Override
    public void openGenome(Genome genome) {
        this.genomeId = genome.getId();
        Arrays.fill(this.roleCounts, 0);
    }

    @Override
    public void recordFeature(Feature feat, String role, int count) {
        int idx = this.getRoleIdx(role);
        if (idx > 0)
            this.roleCounts[idx - 1]++;
    }

    @Override
    public void closeGenome() {
        // Here we write the counts for this genome.
        String counts = Arrays.stream(this.roleCounts).mapToObj(x -> Integer.toString(x)).collect(Collectors.joining("\t"));
        this.print("%s\t%s", this.genomeId, counts);
    }

    @Override
    public void closeReport() {
    }

}
