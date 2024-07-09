package robot;

import robocode.Robot;
import robocode.ScannedRobotEvent;

public class EscapeRobot extends Robot {

    public void run() {
        while (true) {
            turnRadarRight(360);
        }
    }

    public void onScannedRobot(ScannedRobotEvent event) {
        double enemyBearing = event.getBearingRadians();
        double robotHeading = Math.toRadians(getHeading());
        double enemyDistance = event.getDistance();
        double battlefieldWidth = getBattleFieldWidth();
        double battlefieldHeight = getBattleFieldHeight();
        double gunHeading = Math.toRadians(getGunHeading());
        double gunTurnAngle = robocode.util.Utils.normalRelativeAngle(enemyBearing - gunHeading);

        if (getX() < 50 || getX() > getBattleFieldWidth() - 50 ||
                getY() < 50 || getY() > getBattleFieldHeight() - 50) {
            double middleX = getBattleFieldWidth() / 2;
            double middleY = getBattleFieldHeight() / 2;

            double angleToMiddle = Math.atan2(middleX - getX(), middleY - getY());
            double distanceToMiddle = Math.sqrt(Math.pow(middleX - getX(), 2) + Math.pow(middleY - getY(), 2));

            turnRight(Math.toDegrees(angleToMiddle));
            ahead(distanceToMiddle);

        } else {
            double[] escapeVector = calculateEscapeVector(enemyBearing, robotHeading, enemyDistance,
                    battlefieldWidth, battlefieldHeight);

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
        }
    }

    public void onHitByBullet(HitByBulletEvent event){

    }

    public void onHitWall(HitWallEvent event){

    }

    private double[] calculateEscapeVector(double enemyBearing, double robotHeading, double enemyDistance,
                                           double battlefieldWidth, double battlefieldHeight) {
        double horizontalAngle;
        if (enemyBearing > Math.PI / 2 && enemyBearing < 3 * Math.PI / 2) {
            horizontalAngle = Math.PI - enemyBearing;
        } else {
            horizontalAngle = -enemyBearing;
        }

        double verticalAngle;
        if (enemyBearing > Math.PI) {
            verticalAngle = 2 * Math.PI - enemyBearing;
        } else {
            verticalAngle = Math.PI - enemyBearing;
        }

        double x = -Math.cos(horizontalAngle) * enemyDistance;
        double y = -Math.cos(verticalAngle) * enemyDistance;

        double escapeAngle = robotHeading + Math.atan2(x, y);
        double escapeDistance = Math.sqrt(x * x + y * y);

        double xAdjusted = Math.sin(escapeAngle) * escapeDistance;
        double yAdjusted = Math.cos(escapeAngle) * escapeDistance;

        double xWall = Math.min(Math.abs(battlefieldWidth - getX()), getX());
        double yWall = Math.min(Math.abs(battlefieldHeight - getY()), getY());

        xAdjusted += xWall;
        yAdjusted += yWall;

        double xFinal = Math.min(Math.max(50, getX() + xAdjusted), battlefieldWidth - 50);
        double yFinal = Math.min(Math.max(50, getY() + yAdjusted), battlefieldHeight - 50);

        double angleToFinal = Math.atan2(xFinal - getX(), yFinal - getY());

        return new double[]{escapeDistance, Math.toDegrees(angleToFinal)};
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
