/**
 *
 */
package org.theseed.genome.compare;

import java.security.NoSuchAlgorithmException;

import org.theseed.genome.Feature;
import org.theseed.proteins.Function;
import org.theseed.proteins.FunctionMap;

/**
 * This comparison processor counts functional differences in genomes.  Differences in gene-calling are ignored.
 * The intent is to track improvements in the functional assignments.
 *
 * @author Bruce Parrello
 *
 */
public class CompareGenomes extends CompareORFs {

    // FIELDS
    /** good counts */
    int good;
    /** bad counts */
    int bad;
    /** function definition map */
    private FunctionMap funMap;

    /**
     * Initialize the comparison processor.
     *
     * @throws NoSuchAlgorithmException
     */
    public CompareGenomes() throws NoSuchAlgorithmException {
        super();
        this.funMap = new FunctionMap();
    }

    @Override
    protected void both(Feature oldFeature, Feature newFeature) {
        // Compare the functions.
        String oldFun = oldFeature.getPegFunction();
        String newFun = newFeature.getPegFunction();
        Function fun = funMap.findOrInsert(oldFun);
        Function other = funMap.getByName(newFun);
        // Count good if the functions are equal.
        if (other != null && other.equals(fun))
            good++;
        else
            bad++;
    }

    @Override
    protected void newOnly(Feature newFeature) {
    }

    @Override
    protected void oldOnly(Feature oldFeature) {
    }

    @Override
    protected void initCompareData() {
        this.good = 0;
        this.bad = 0;
    }

    /**
     * @return the percent of features that were good
     */
    public double percent() {
        double retVal = 0.0;
        if (this.good > 0) {
            retVal = good * 100.0 / (good + bad);
        }
        return retVal;
    }

    /**
     * @return the nubmer of good matches
     */
    public int getGood() {
        return this.good;
    }

    /**
     * @return the number of bad matches
     */
    public int getBad() {
        return this.bad;
    }

}
