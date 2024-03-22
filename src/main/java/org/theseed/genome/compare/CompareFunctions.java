/**
 *
 */
package org.theseed.genome.compare;

import java.security.NoSuchAlgorithmException;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import org.theseed.counters.CountMap;
import org.theseed.genome.Feature;
import org.theseed.proteins.Function;
import org.theseed.proteins.FunctionMap;
import org.theseed.stats.QualityCountMap;

/**
 * This class tracks comparisons between functional annotations.  We expect most functions to map to themselves, but we
 * track the ones that map to something else.  Only features that are in both the new genome AND the old genome are
 * tracked.  For those, we map the old annotation to the new one.  The mapping will be by function ID, and we will
 * use a simple counter to track the identity mappings, and a CountMap for the non-identity mappings.
 *
 * @author Bruce Parrello
 *
 */
public class CompareFunctions extends CompareORFs {

    // FIELDS
    /** function ID map */
    private FunctionMap funMap;
    /** matching-annotation (good) and missed-annotation (bad) counts */
    private QualityCountMap<String> matchCounts;
    /** map of wrong-annotation counts */
    private Map<String, CountMap<String>> missCountMap;
    /** comparator for sorting results */
    private Comparator<Function> comparator;

    /**
     * Initialize the maps.
     *
     * @throws NoSuchAlgorithmException
     */
    public CompareFunctions() throws NoSuchAlgorithmException {
        super();
        this.funMap = new FunctionMap();
        this.matchCounts = new QualityCountMap<String>();
        this.missCountMap = new HashMap<String, CountMap<String>>(5000);
        this.comparator = this.new FunctionCompare();
    }

    @Override
    protected void both(Feature oldFeature, Feature newFeature) {
        // Get the function IDs for the old and new features.
        String oldFun = this.getFunctionId(oldFeature.getFunction());
        String newFun = this.getFunctionId(newFeature.getFunction());
        if (oldFun.contentEquals(newFun))
            this.matchCounts.setGood(oldFun);
        else {
            // Here a function has changed its name.  Denote it is a possible replacement for the old function.
            CountMap<String> missCounts = this.missCountMap.computeIfAbsent(oldFun, f -> new CountMap<String>());
            missCounts.count(newFun);
            this.matchCounts.setBad(oldFun);
        }
    }

    /**
     * @return the ID corresponding to the specified functional assignment
     *
     * @param function	functional assignment whose ID is desired
     */
    private String getFunctionId(String function) {
        Function fun = this.funMap.findOrInsert(function);
        return fun.getId();
    }

    /**
     * @return the count map for the specified function ID, or NULL if there were no misses
     *
     * @param funId		ID of the function of interest
     */
    public CountMap<String> getMissCounts(String funId) {
        return this.missCountMap.get(funId);
    }

    /**
     * @return the match count for the specified function ID
     *
     * @param funId		ID of the function of interest
     */
    public int getMatchCount(String funId) {
        return this.matchCounts.good(funId);
    }

    /**
     * @return the total count for the specified function ID
     *
     * @param funId		ID of the function of interest
     */
    public int getTotalCount(String funId) {
        return this.matchCounts.good(funId) + this.matchCounts.bad(funId);
    }

    /**
     * This is a comparator for sorting functions by lowest to highest good count.
     */
    private class FunctionCompare implements Comparator<Function> {

        @Override
        public int compare(Function o1, Function o2) {
            int retVal = CompareFunctions.this.matchCounts.good(o1.getId()) -
                    CompareFunctions.this.matchCounts.good(o2.getId());
            if (retVal == 0)
                retVal = o1.getName().compareTo(o2.getName());
            return retVal;
        }

    }

    /**
     * @return a collection of the functions with miss counts
     */
    public SortedSet<Function> getMissFunctions() {
        SortedSet<Function> retVal = new TreeSet<Function>(this.comparator);
        this.missCountMap.keySet().stream().map(f -> this.funMap.getItem(f)).forEach(x -> retVal.add(x));
        return retVal;
    }

    @Override
    protected void newOnly(Feature newFeature) {
    }

    @Override
    protected void oldOnly(Feature oldFeature) {
    }

    @Override
    protected void initCompareData() {
    }

    /**
     * @return the name of the function with the specified ID
     *
     * @param funId		function of interest
     * @return
     */
    public String getName(String funId) {
        return this.funMap.getName(funId);
    }

}
