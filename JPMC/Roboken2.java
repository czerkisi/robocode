package JPMC;

import robocode.Robot;
import robocode.HitRobotEvent;
import robocode.ScannedRobotEvent;

import java.awt.*;

/**
 * SpinBot - a sample robot by Mathew Nelson.
 * <p>
 * Moves in a circle, firing hard when an enemy is detected.
 *
 * @auther Mathew A. Nelson (original)
 * @auther Flemming N. Larsen (contributor)
 */
public class Roboken2 extends Robot {

    /**
     * SpinBot's run method - Circle
     */
    public void run() {
        // Set colors
        setColors(Color.blue, Color.blue, Color.black, Color.yellow, Color.white);

        // Loop forever
        while (true) {
            // Turn the robot and move ahead indefinitely
            turnRight(10000);
            ahead(10000);
        }
    }

    /**
     * onScannedRobot: Fire hard!
     */
    public void onScannedRobot(ScannedRobotEvent e) {
        fire(3);
    }

    /**
     * onHitRobot:  If it's our fault, we'll stop turning and moving,
     * so we need to turn again to keep spinning.
     */
    public void onHitRobot(HitRobotEvent e) {
        if (e.getBearing() > -10 && e.getBearing() < 10) {
            fire(3);
        }
        if (e.isMyFault()) {
            turnRight(10);
        }
    }
}