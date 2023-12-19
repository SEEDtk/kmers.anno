/**
 *
 */
package org.theseed.reports;

import java.util.Arrays;
import java.util.List;

import org.theseed.genome.Feature;
import org.theseed.genome.SubsystemRow;

/**
 * This report outputs every feature pair, and includes both the annotations and the subsystem
 * data.
 */
public class FullCompareAnnotationReporter extends AnnotationReporter {

    // FIELDS
    /** output buffer */
    private final Object[] fields;

    /**
     * Construct a new full-comparison reporter.
     */
    public FullCompareAnnotationReporter() {
        this.fields = new Object[11];
    }

    @Override
    protected void start(IParms processor) {
        this.writeHeader("fid", "old_annotation", "old_subsystem", "old_subclass1",
                "old_subclass2", "old_subclass3", "new_annotation", "new_subsystem",
                "new_subclass1", "new_subclass2", "new_subclass3");
    }

    @Override
    public void processFeature(Feature oldFeat, Feature newFeat) {
        String fid = oldFeat.getId();
        String oldAnno = oldFeat.getPegFunction();
        String newAnno = newFeat.getPegFunction();
        // Now we need to compare the subsystems.  We get a list of associated subsystems for each version.
        var oldSubs = oldFeat.getSubsystemRows();
        var newSubs = newFeat.getSubsystemRows();
        if (oldSubs.isEmpty() && newSubs.isEmpty()) {
            // Here there are no subsystems.  Show the annotations only.
            this.writeData(fid, oldAnno, null, null, null, null, newAnno, null, null, null, null);
        } else {
            var oldIter = oldSubs.iterator();
            var newIter = newSubs.iterator();
            while (oldIter.hasNext() && newIter.hasNext()) {
                // Clear out the old data line.
                Arrays.fill(this.fields, null);
                // Start filling in the new data.
                this.fields[0] = fid;
                this.fields[1] = oldAnno;
                if (oldIter.hasNext()) {
                    SubsystemRow oldRow = oldIter.next();
                    fillSubData(oldRow, this.fields, 2);
                }
                this.fields[6] = newAnno;
                if (newIter.hasNext()) {
                    SubsystemRow newRow = newIter.next();
                    fillSubData(newRow, this.fields, 7);
                }
                // Write the line.
                this.writeData(this.fields);
            }
        }
    }

    /**
     * Fill in the subsystem name and classes.
     *
     * @param subRow	subsystem data
     * @param fields	output field buffer
     * @param i			position for subsystem name
     */
    protected static void fillSubData(SubsystemRow subRow, Object[] fields, int i) {
        fields[i] = subRow.getName();
        List<String> classes = subRow.getClassifications();
        int n = classes.size();
        if (n > 3) n = 3;
        for (int j = 0; j < n; j++)
            fields[i + j] = classes.get(j);
    }

    @Override
    protected void finish() {
        //  No cleanup needed.
    }

}
