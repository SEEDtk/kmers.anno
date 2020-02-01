/**
 *
 */
package org.theseed.locations;

import java.util.ArrayList;
import java.util.Map;
import org.theseed.locations.FramedLocationLists.Report;

import java.util.HashMap;
import java.util.Iterator;

/**
 * This is a map of location lists, one per target per frame.  The targets are usually feature IDs; however, they
 * can be any string.  The map can be thought of as a mapping from [frame, target] pairs to SortedLocationList
 * objects.
 *
 * The location lists are organized by frame and then target within frame.  The order of organization does not
 * matter to the caller, it is simply for implementation convenience.
 *
 * Methods are provided to insert new locations and to iterate through the lists.
 *
 * @author Bruce Parrello
 *
 */
public class FramedLocationLists implements Iterable<Report> {

    // FIELDS
    /** master array of target/list maps, the array index is frame */
    private ArrayList<Map<String, SortedLocationList>> master;
    /** number of items in all the maps */
    private int count;

    /**
     * This class is used to report on a target containing a location list.
     * The client can iterate through these to get all the lists in the
     * map.
     */
    public class Report {

        // FIELDS
        private String targetId;
        private SortedLocationList list;

        private Report(String target, SortedLocationList list) {
            this.targetId = target;
            this.list = list;
        }

        /**
         * @return the ID of the target relevant to this location list
         */
        public String getId() {
            return targetId;
        }

        /**
         * @return the list of locations for the target
         */
        public SortedLocationList getList() {
            return list;
        }

    }

    /**
     * Iterator for this class
     */
    public class FramedLocationListsIterator implements Iterator<Report> {

        // FIELDS
        /** iterator through the array list */
        private Iterator<Map<String, SortedLocationList>> frameIter;
        /** iterator through the keys of the current map */
        private Iterator<String> mapIter;
        /** current map of targets to location lists */
        private Map<String, SortedLocationList> currentMap;

        private FramedLocationListsIterator() {
            this.frameIter = FramedLocationLists.this.master.iterator();
            this.currentMap = this.frameIter.next();
            this.mapIter = currentMap.keySet().iterator();
        }

        @Override
        public boolean hasNext() {
            // Insure we are on a nonempty map.
            this.fixFrame();
            return this.mapIter.hasNext();
        }

        @Override
        public Report next() {
            // Advance the frame if the current map is done.
            this.fixFrame();
            // Get the next target ID.
            String target = this.mapIter.next();
            // Return a report containing the target ID and the location list.
            return new Report(target, this.currentMap.get(target));
        }

        /**
         * Advance the frame until we are at the end or on a nonempty map.  Note
         * that the last map is ALWAYS empty. (This is an artifact.)
         */
        private void fixFrame() {
            while (this.frameIter.hasNext() && ! this.mapIter.hasNext()) {
                this.currentMap = this.frameIter.next();
                this.mapIter = this.currentMap.keySet().iterator();
            }
        }

    }

    /**
     * Construct an empty mapping.
     */
    public FramedLocationLists() {
        // Create the array.
        this.master = new ArrayList<Map<String, SortedLocationList>>(Frame.nFrames);
        // Fill in blank maps for all the frames.
        for (Frame frm : Frame.all) {
            this.master.add(frm.getIdx(), new HashMap<String, SortedLocationList>());
        }
        this.count = 0;
    }

    @Override
    public Iterator<Report> iterator() {
        return this.new FramedLocationListsIterator();
    }

    /**
     * Clear all the entries in the object.
     */
    public void clear() {
        for (Map<String, SortedLocationList> map : this.master) {
            map.clear();
        }
        this.count = 0;
    }

    /**
     * @return the total number of connections
     */
    public int size() {
        return count;
    }

    /**
     * Connect a location to a specific target.
     *
     * @param target	ID of the target
     * @param loc		location to insert
     */
    public void connect(String target, Location loc) {
        // Compute the location's frame.
        Frame locFrame = loc.getFrame();
        // Get the map for that frame.
        Map<String, SortedLocationList> map = this.master.get(locFrame.getIdx());
        // Insure there is a location list for the target.
        SortedLocationList locList = map.get(target);
        if (locList == null) {
            locList = new SortedLocationList();
            map.put(target, locList);
        }
        // Add this location into the list.
        locList.add(loc);
        // Increment the count.
        this.count++;
    }

}
