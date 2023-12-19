/**
 *
 */
package org.theseed.reports;

import java.util.Arrays;

import org.theseed.genome.Feature;
import org.theseed.genome.SubsystemRow;


/**
 * This report only outputs features for which the original role was hypothetical.
 */
public class NewRoleAnnotationReporter extends AnnotationReporter {

    // FIELDS
    /** output buffer */
    private final Object[] fields;

    /**
     * Construct a new full-comparison reporter.
     */
    public NewRoleAnnotationReporter() {
        this.fields = new Object[7];
    }

    @Override
    protected void start(IParms processor) {
        this.writeHeader("fid", "old_annotation", "new_annotation", "new_subsystem",
                "new_subclass1", "new_subclass2", "new_subclass3");
    }

    @Override
    public void processFeature(Feature oldFeat, Feature newFeat) {
        String fid = oldFeat.getId();
        String oldAnno = oldFeat.getPegFunction();
        String newAnno = newFeat.getPegFunction();
        if (oldAnno.equals("hypothetical protein") && ! oldAnno.equals(newAnno)) {
            // Now we need to output this annotation.  For the new one, we need subsystem data.
            var newSubs = newFeat.getSubsystemRows();
            if (newSubs.isEmpty()) {
                // Here there are no subsystems.  Show the annotations only.
                this.writeData(fid, oldAnno, newAnno, null, null, null, null);
            } else {
                var newIter = newSubs.iterator();
                while (newIter.hasNext()) {
                    // Clear out the old data line.
                    Arrays.fill(this.fields, null);
                    // Start filling in the new data.
                    this.fields[0] = fid;
                    this.fields[1] = oldAnno;
                    this.fields[2] = newAnno;
                    SubsystemRow newRow = newIter.next();
                    FullCompareAnnotationReporter.fillSubData(newRow, this.fields, 3);
                    // Write the line.
                    this.writeData(this.fields);
                }
            }
        }
    }

    @Override
    protected void finish() {
        // No resources to release.
    }

}
