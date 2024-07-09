package JPMC;

import robocode.Robot;
import robocode.ScannedRobotEvent;

public class Roboken2 extends Robot {

    public void run() {
        while (true) {
            turnRadarRight(360);  // Continuous radar scan
        }
    }

    public void onScannedRobot(ScannedRobotEvent event) {
        double enemyBearing = event.getBearingRadians();
        double enemyDistance = event.getDistance();
        double gunTurnAngle = robocode.util.Utils.normalRelativeAngle(enemyBearing - Math.toRadians(getGunHeading()));

        if (isNearWall()) {
            moveToCenter();
        } else {
            double[] escapeVector = calculateEscapeVector(enemyBearing, enemyDistance);
            double distanceToFinal = escapeVector[0];
            double turnAngle = escapeVector[1];

            double wallWeight = calculateWallWeight();
            double obstacleWeight = calculateObstacleWeight(enemyDistance);

            ahead(distanceToFinal * (1 - wallWeight) * (1 - obstacleWeight));
            turnRight(turnAngle);
        }

        if (enemyDistance < 100) {
            turnGunRight(Math.toDegrees(gunTurnAngle));
            fire(3);
        } else {
            // Predictive targeting
            double bulletPower = Math.min(3.0, getEnergy() / 10);
            double bulletSpeed = 20 - 3 * bulletPower;
            double time = enemyDistance / bulletSpeed;
            double futureX = event.getVelocity() * time * Math.sin(event.getHeadingRadians());
            double futureY = event.getVelocity() * time * Math.cos(event.getHeadingRadians());
            double futureBearing = Math.atan2(futureX, futureY);

            turnGunRight(Math.toDegrees(robocode.util.Utils.normalRelativeAngle(futureBearing - Math.toRadians(getGunHeading()))));
            fire(bulletPower);
        }
    }

    private boolean isNearWall() {
        return getX() < 50 || getX() > getBattleFieldWidth() - 50 ||
                getY() < 50 || getY() > getBattleFieldHeight() - 50;
    }

    private void moveToCenter() {
        double middleX = getBattleFieldWidth() / 2;
        double middleY = getBattleFieldHeight() / 2;
        double angleToMiddle = Math.atan2(middleX - getX(), middleY - getY());
        double distanceToMiddle = Math.sqrt(Math.pow(middleX - getX(), 2) + Math.pow(middleY - getY(), 2));

        turnRight(Math.toDegrees(angleToMiddle) - getHeading());
        ahead(distanceToMiddle);
    }

    private double[] calculateEscapeVector(double enemyBearing, double enemyDistance) {
        double escapeAngle = Math.toRadians(getHeading()) + enemyBearing + Math.PI / 2;
        double escapeDistance = 100;  // Arbitrary escape distance

        return new double[]{escapeDistance, Math.toDegrees(escapeAngle)};
    }

    private double calculateWallWeight() {
        double wallWeight = 0;
        double xDistanceToWall = Math.min(getX(), getBattleFieldWidth() - getX());
        double yDistanceToWall = Math.min(getY(), getBattleFieldHeight() - getY());

        if (xDistanceToWall < 100) {
            wallWeight += (100 - xDistanceToWall) / 100;
        }
        if (yDistanceToWall < 100) {
            wallWeight += (100 - yDistanceToWall) / 100;
        }

        return wallWeight;
    }

    private double calculateObstacleWeight(double enemyDistance) {
        double obstacleWeight = 0;

        if (enemyDistance < 200) {
            obstacleWeight = (200 - enemyDistance) / 200;
        }

        return obstacleWeight;
    }
}
