/**
 *
 */
package org.theseed.proteins.kmers;

import junit.framework.TestCase;
import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;
/**
 * @author Bruce Parrello
 *
 */
public class RoleTests extends TestCase {

    public void testRoleCounter() {
        RoleCounter counter = new RoleCounter("roleA");
        assertThat(counter.getRoleId(), equalTo("roleA"));
        assertThat(counter.getGoodCount(), equalTo(0));
        assertThat(counter.getBadCount(), equalTo(0));
        assertThat(counter.isGood(), equalTo(true));
        counter.count("roleA");
        assertThat(counter.getRoleId(), equalTo("roleA"));
        assertThat(counter.getGoodCount(), equalTo(1));
        assertThat(counter.getBadCount(), equalTo(0));
        assertThat(counter.isGood(), equalTo(true));
        counter.count("roleA");
        assertThat(counter.getRoleId(), equalTo("roleA"));
        assertThat(counter.getGoodCount(), equalTo(2));
        assertThat(counter.getBadCount(), equalTo(0));
        assertThat(counter.isGood(), equalTo(true));
        counter.count("roleB");
        assertThat(counter.getRoleId(), equalTo("roleA"));
        assertThat(counter.getGoodCount(), equalTo(2));
        assertThat(counter.getBadCount(), equalTo(1));
        assertThat(counter.isGood(), equalTo(false));
    }
}
