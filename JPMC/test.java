public class EscapeRobot extends AdvancedRobot {

    // Run method
    public void run() {
        setAdjustGunForRobotTurn(true);
        setAdjustRadarForGunTurn(true);
        setAdjustRadarForRobotTurn(true);
        while (true) {
            turnRadarRight(360);  // Continuous radar scan
        }
    }

    // When the robot scans another robot
    public void onScannedRobot(ScannedRobotEvent event) {
        double enemyBearing = event.getBearingRadians();
        double gunTurnAngle = robocode.util.Utils.normalRelativeAngle(enemyBearing - getGunHeadingRadians());

        // Improved predictive targeting
        double enemyVelocity = event.getVelocity();
        double enemyHeading = event.getHeadingRadians();
        double bulletPower = Math.min(3.0, getEnergy());

        double predictedX = getX() + event.getDistance() * Math.sin(enemyBearing);
        double predictedY = getY() + event.getDistance() * Math.cos(enemyBearing);

        // Simple linear prediction
        double predictionTime = event.getDistance() / (20 - 3 * bulletPower);
        predictedX += Math.sin(enemyHeading) * enemyVelocity * predictionTime;
        predictedY += Math.cos(enemyHeading) * enemyVelocity * predictionTime;

        double aimAngle = robocode.util.Utils.normalRelativeAngle(Math.atan2(predictedX - getX(), predictedY - getY()));
        setTurnGunRightRadians(aimAngle - getGunHeadingRadians());
        if (getGunHeat() == 0) {
            setFire(bulletPower);
        }

        // Improved movement with non-linear pattern
        setAhead(100 * (Math.random() - 0.5));
        setTurnRightRadians(robocode.util.Utils.normalRelativeAngle(Math.PI / 2 - enemyBearing));

        // Execute all pending actions
        execute();
    }
}
