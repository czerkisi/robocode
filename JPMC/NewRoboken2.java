public class SpinBot extends Robot {

    @Override
    public void run() {
        // Set colors
        setBodyColor(Color.blue);
        setGunColor(Color.blue);
        setRadarColor(Color.black);
        setScanColor(Color.yellow);

        // Get the dimensions of the battlefield
        double battlefieldWidth = getBattleFieldWidth();
        double battlefieldHeight = getBattleFieldHeight();

        // Calculate the center point of the battlefield
        double centerX = battlefieldWidth / 2;
        double centerY = battlefieldHeight / 2;

        // Calculate the radius for the figure-8 pattern
        double radius = Math.min(battlefieldWidth, battlefieldHeight) / 4;

        // Loop forever
        while (true) {
            moveInFigure8(centerX, centerY, radius);
        }
    }

    private void moveInFigure8(double centerX, double centerY, double radius) {
        // Move in a circle to the right
        for (int i = 0; i < 180; i++) {
            turnRight(2);
            ahead(Math.PI * radius / 180);
        }
        // Move in a circle to the left
        for (int i = 0; i < 180; i++) {
            turnLeft(2);
            ahead(Math.PI * radius / 180);
        }
    }

    @Override
    public void onScannedRobot(ScannedRobotEvent e) {
        fire(3);
    }

    @Override
    public void onHitRobot(HitRobotEvent e) {
        if (e.getBearing() > -10 && e.getBearing() < 10) {
            fire(3);
        }
        if (e.isMyFault()) {
            turnRight(10);
        }
    }
}
