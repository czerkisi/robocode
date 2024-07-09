package JPMC;

import robocode.Robot;

public class Roboken2 extends Robot {

    public void run() {
        // Loop to continuously perform the figure-eight movement
        while (true) {
            moveInFigureEight();
        }
    }

    private void moveInFigureEight() {
        // Move in one loop of the figure-eight
        for (int i = 0; i < 2; i++) {
            // Move in a circular pattern
            for (int j = 0; j < 36; j++) {
                ahead(10);
                turnRight(10);
            }
            // Reverse the direction for the second loop
            for (int j = 0; j < 36; j++) {
                ahead(10);
                turnLeft(10);
            }
        }
    }
}
