import robocode.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class WinningRobot extends Robot {
    private static final int MAX_COORDINATES = 100;
    private static final int MIN_COORDINATES = 20;
    private static final double EPSILON = 1e-6;

    private Map<String, List<Coordinate>> robotMovements;

    @Override
    public void run() {
        robotMovements = new HashMap<>();

        while (true) {
            avoidRobots();
            avoidWalls();
            turnRadarRight(360);
        }
    }

    private void avoidRobots() {
        double safeDistance = 150; // Minimum safe distance from other robots

        // Find the closest robot
        double closestDistance = Double.MAX_VALUE;
        double closestBearing = 0;

        for (String enemyName : robotMovements.keySet()) {
            List<Coordinate> coordinates = robotMovements.get(enemyName);
            if (coordinates.isEmpty()) {
                continue;
            }

            Coordinate latestCoordinate = coordinates.get(coordinates.size() - 1);
            double distance = getDistanceBetween(new Coordinate(getX(), getY()), latestCoordinate);
            if (distance < closestDistance) {
                closestDistance = distance;
                double deltaX = latestCoordinate.getX() - getX();
                double deltaY = latestCoordinate.getY() - getY();
                closestBearing = Math.atan2(deltaX, deltaY);
            }
        }

        // Move away from the closest robot if within safe distance
        if (closestDistance < safeDistance) {
            double angleToRobot = Math.toDegrees(closestBearing) - getHeading();
            turnRight(angleToRobot);
            ahead(safeDistance - closestDistance);
        }
    }

    private void avoidWalls() {
        double wallMargin = 50; // Margin to consider near the wall

        // Check if near a wall
        if (isNearWall()) {
            double wallBearing = getBearingToWall();
            double wallTurn = robocode.util.Utils.normalRelativeAngleDegrees(wallBearing - getHeading());

            if (wallTurn >= 0) {
                turnRight(wallTurn + 90);
            } else {
                turnLeft(Math.abs(wallTurn) + 90);
            }

            ahead(100);
        }
    }

    private boolean isNearWall() {
        double wallMargin = 50; // Margin to consider near the wall
        return (getX() <= wallMargin || getY() <= wallMargin || getBattleFieldWidth() - getX() <= wallMargin
                || getBattleFieldHeight() - getY() <= wallMargin);
    }

    private double getBearingToWall() {
        double x = getX();
        double y = getY();
        double battlefieldWidth = getBattleFieldWidth();
        double battlefieldHeight = getBattleFieldHeight();

        double northBearing = Math.toDegrees(Math.atan2(battlefieldHeight - y, x));
        double southBearing = Math.toDegrees(Math.atan2(-y, x));
        double eastBearing = Math.toDegrees(Math.atan2(battlefieldWidth - x, y));
        double westBearing = Math.toDegrees(Math.atan2(x, y));

        double wallBearing = Math.min(Math.min(northBearing, southBearing), Math.min(eastBearing, westBearing));
        return wallBearing;
    }

    private double calculateTimeToHit(double targetX, double targetY) {
        double bulletPower = 2;  // Set the desired firepower of the bullet

        // Calculate the distance between the robot and the target coordinates
        double deltaX = targetX - getX();
        double deltaY = targetY - getY();
        double distance = Math.sqrt(deltaX * deltaX + deltaY * deltaY);

        // Calculate the angle between the robot's gun heading and the target
        double absoluteBearing = Math.atan2(deltaX, deltaY);
        double gunHeading = Math.toRadians(getGunHeading());
        double gunTurnRate = Math.toRadians(20);  // Maximum rate of gun rotation in radians/turn

        double deltaAngle = robocode.util.Utils.normalRelativeAngle(absoluteBearing - gunHeading);
        double timeToTurn = Math.abs(deltaAngle) / gunTurnRate;

        // Calculate the time it takes for the bullet to reach the target
        double bulletVelocity = 20 - 3 * bulletPower;
        double timeToHit = distance / bulletVelocity + timeToTurn;

        return timeToHit;
    }

    private void aimAndShoot(double bulletPower, RegressionAnalyzer analyzer) {
        // Get the predicted future coordinates of the enemy robot
        Coordinate currentCoordinate = analyzer.predict(0);
        double predictionTime = calculateTimeToHit(currentCoordinate.getX(), currentCoordinate.getY());
        Coordinate predictedCoordinates = analyzer.predict((int) predictionTime);

        // Calculate the angle between the robot's gun heading and the predicted coordinates
        double deltaX = predictedCoordinates.getX() - getX();
        double deltaY = predictedCoordinates.getY() - getY();
        double absoluteBearing = Math.atan2(deltaX, deltaY);
        double gunHeading = Math.toRadians(getGunHeading());
        double gunTurnRate = Math.toRadians(20);  // Maximum rate of gun rotation in radians/turn

        double deltaAngle = robocode.util.Utils.normalRelativeAngle(absoluteBearing - gunHeading);
        double timeToTurn = Math.abs(deltaAngle) / gunTurnRate;

        // Adjust the prediction time and gun angle iteratively until the bullet hits the enemy
        double bulletVelocity = 20 - 3 * bulletPower;
        double timeToHit = predictionTime;  // Initial prediction time

        while (timeToHit > predictionTime) {
            double distance = getDistanceBetween(predictedCoordinates, new Coordinate(getX(), getY()));
            timeToHit = distance / bulletVelocity + timeToTurn;

            // Check if the gun can be turned and fired in time
            if (timeToHit <= predictionTime) {
                // Turn the gun to the predicted coordinates
                turnGunRight(Math.toRadians(deltaAngle));

                // Fire the gun
                fire(bulletPower);
                break;
            }

            // Adjust the prediction time and recalculate the predicted coordinates
            predictionTime += 1;  // Increase the prediction time by 1 tick
            predictedCoordinates = analyzer.predict((int) predictionTime);

            // Recalculate the angle between the robot's gun heading and the new predicted coordinates
            deltaX = predictedCoordinates.getX() - getX();
            deltaY = predictedCoordinates.getY() - getY();
            absoluteBearing = Math.atan2(deltaX, deltaY);
            deltaAngle = robocode.util.Utils.normalRelativeAngle(absoluteBearing - gunHeading);
            timeToTurn = Math.abs(deltaAngle) / gunTurnRate;
        }
    }



    @Override
    public void onScannedRobot(ScannedRobotEvent event) {
        double enemyDistance = event.getDistance();
        double enemyHeading = getHeading() + event.getBearingRadians();

        // Calculate the coordinates of the scanned robot
        double enemyX = getX() + Math.sin(enemyHeading) * enemyDistance;
        double enemyY = getY() + Math.cos(enemyHeading) * enemyDistance;

        String enemyName = event.getName();

        // Update or create the list of coordinates for the scanned robot
        List<Coordinate> coordinateList = robotMovements.get(enemyName);
        if (coordinateList == null) {
            coordinateList = new ArrayList<>();
            robotMovements.put(enemyName, coordinateList);
        }

        // Store the current coordinates of the scanned robot
        coordinateList.add(new Coordinate(enemyX, enemyY));

        // Trim the list to a maximum of 100 coordinates
        if (coordinateList.size() > MAX_COORDINATES) {
            coordinateList.remove(0);
        }
        if (coordinateList.size() >= MIN_COORDINATES){
            determineAction(enemyName);
        }
    }

    private double getDistanceBetween(Coordinate coord1, Coordinate coord2) {
        double deltaX = coord2.getX() - coord1.getX();
        double deltaY = coord2.getY() - coord1.getY();

        return Math.sqrt(deltaX * deltaX + deltaY * deltaY);
    }

    private void determineAction(String enemyName) {
        List<Coordinate> coordinates = robotMovements.get(enemyName);
        RegressionAnalyzer analyzer = new RegressionAnalyzer(coordinates);
        Coordinate currentCoordinate = coordinates.get(coordinates.size() - 1);
        if (getDistanceBetween(currentCoordinate, new Coordinate(getX(), getY())) < 20){
            aimAndShoot(5, analyzer);
        }
    }

    private static class Coordinate {
        private double x;
        private double y;

        public Coordinate(double x, double y) {
            this.x = x;
            this.y = y;
        }

        public double getX() {
            return x;
        }

        public double getY() {
            return y;
        }
    }

    public static class RegressionAnalyzer {
        private List<Coordinate> coordinates;
        private RegressionModel bestFitModel;

        public RegressionAnalyzer(List<Coordinate> coordinates) {
            this.coordinates = coordinates;
            this.bestFitModel = determineBestFitModel();
        }

        public Coordinate predict(int n) {
            return bestFitModel.predict(n);
        }

        public boolean hasPattern() {
            return bestFitModel != null;
        }

        private RegressionModel determineBestFitModel() {
            RegressionModel bestModel = null;
            double bestModelError = Double.MAX_VALUE;

            RegressionModel[] models = {new LinearRegressionModel(), new QuadraticRegressionModel(), new CircularRegressionModel(), new SpiralRegressionModel(), new StillRegressionModel()};
            for (RegressionModel model : models) {
                double modelError = calculateModelError(model);

                if (modelError < bestModelError) {
                    bestModelError = modelError;
                    bestModel = model;
                }
            }

            return bestModel;
        }

        private double calculateModelError(RegressionModel model) {
            double error = 0.0;

            for (int i = 0; i < coordinates.size(); i++) {
                Coordinate coordinate = coordinates.get(i);
                Coordinate predictedCoordinate = model.predict(i);
                error += model.calculateError(predictedCoordinate, coordinate);
            }

            return error;
        }

        public abstract class RegressionModel {
            protected List<Coordinate> coordinates;

            public abstract void train(List<Coordinate> coordinates);
            public abstract Coordinate predict(int n);

            protected double calculateError(Coordinate predictedCoordinate, Coordinate actualCoordinate) {
                double xDiff = actualCoordinate.getX() - predictedCoordinate.getX();
                double yDiff = actualCoordinate.getY() - predictedCoordinate.getY();

                return Math.sqrt(xDiff * xDiff + yDiff * yDiff);
            }
        }

        public class LinearRegressionModel extends RegressionModel {
            private double slope;
            private double intercept;

            @Override
            public void train(List<Coordinate> coordinates) {
                this.coordinates = coordinates;
                calculateParameters();
            }

            @Override
            public Coordinate predict(int n) {
                if (slope == 0) {
                    throw new IllegalStateException("Model has not been trained yet.");
                }

                double x = n;
                double y = slope * x + intercept;

                return new Coordinate(x, y);
            }

            private void calculateParameters() {
                int numPoints = coordinates.size();

                double sumN = 0;
                double sumX = 0;
                double sumY = 0;
                double sumNX = 0;
                double sumNX2 = 0;

                for (Coordinate coordinate : coordinates) {
                    double n = coordinates.indexOf(coordinate);
                    double x = coordinate.getX();
                    double y = coordinate.getY();

                    sumN += n;
                    sumX += x;
                    sumY += y;
                    sumNX += n * x;
                    sumNX2 += n * n * x;
                }

                double numerator = (numPoints * sumNX) - (sumN * sumX);
                double denominator = (numPoints * sumNX2) - (sumN * sumN);

                slope = numerator / denominator;
                intercept = (sumY - slope * sumX) / numPoints;
            }
        }

        public class QuadraticRegressionModel extends RegressionModel {
            private double a;
            private double b;
            private double c;

            @Override
            public void train(List<Coordinate> coordinates) {
                this.coordinates = coordinates;
                calculateParameters();
            }

            @Override
            public Coordinate predict(int n) {
                if (a == 0 && b == 0 && c == 0) {
                    throw new IllegalStateException("Model has not been trained yet.");
                }

                double x = n;
                double y = a * Math.pow(x, 2) + b * x + c;

                return new Coordinate(x, y);
            }

            private void calculateParameters() {
                int numPoints = coordinates.size();

                double sumX = 0;
                double sumX2 = 0;
                double sumX3 = 0;
                double sumX4 = 0;
                double sumY = 0;
                double sumXY = 0;
                double sumX2Y = 0;

                for (Coordinate coordinate : coordinates) {
                    double x = coordinate.getX();
                    double y = coordinate.getY();

                    sumX += x;
                    sumX2 += Math.pow(x, 2);
                    sumX3 += Math.pow(x, 3);
                    sumX4 += Math.pow(x, 4);
                    sumY += y;
                    sumXY += x * y;
                    sumX2Y += Math.pow(x, 2) * y;
                }

                double[][] matrix = {
                        { numPoints, sumX, sumX2 },
                        { sumX, sumX2, sumX3 },
                        { sumX2, sumX3, sumX4 }
                };

                double[] constants = { sumY, sumXY, sumX2Y };

                GaussianEliminationSolver solver = new GaussianEliminationSolver(matrix, constants);
                double[] coefficients = solver.solve();

                a = coefficients[0];
                b = coefficients[1];
                c = coefficients[2];
            }

            public class GaussianEliminationSolver {
                private double[][] matrix;
                private double[] constants;

                public GaussianEliminationSolver(double[][] matrix, double[] constants) {
                    this.matrix = matrix;
                    this.constants = constants;
                }

                public double[] solve() {
                    int numRows = matrix.length;
                    int numCols = matrix[0].length;

                    // Forward elimination
                    for (int pivot = 0; pivot < numRows - 1; pivot++) {
                        for (int row = pivot + 1; row < numRows; row++) {
                            double factor = matrix[row][pivot] / matrix[pivot][pivot];
                            constants[row] -= factor * constants[pivot];
                            for (int col = pivot; col < numCols; col++) {
                                matrix[row][col] -= factor * matrix[pivot][col];
                            }
                        }
                    }

                    // Back substitution
                    double[] solution = new double[numRows];
                    for (int row = numRows - 1; row >= 0; row--) {
                        double sum = 0;
                        for (int col = row + 1; col < numCols - 1; col++) {
                            sum += matrix[row][col] * solution[col];
                        }
                        solution[row] = (constants[row] - sum) / matrix[row][row];
                    }

                    return solution;
                }
            }
        }
        public class CircularRegressionModel extends RegressionModel {
            private double centerX;
            private double centerY;
            private double radius;

            @Override
            public void train(List<Coordinate> coordinates) {
                this.coordinates = coordinates;
                calculateParameters();
            }

            @Override
            public Coordinate predict(int n) {
                if (centerX == 0 && centerY == 0 && radius == 0) {
                    throw new IllegalStateException("Model has not been trained yet.");
                }

                double theta = Math.toRadians(n);
                double x = centerX + radius * Math.cos(theta);
                double y = centerY + radius * Math.sin(theta);

                return new Coordinate(x, y);
            }

            private void calculateParameters() {
                int numPoints = coordinates.size();

                double sumX = 0;
                double sumY = 0;
                double sumX2 = 0;
                double sumY2 = 0;
                double sumXY = 0;

                for (Coordinate coordinate : coordinates) {
                    double x = coordinate.getX();
                    double y = coordinate.getY();

                    sumX += x;
                    sumY += y;
                    sumX2 += x * x;
                    sumY2 += y * y;
                    sumXY += x * y;
                }

                double centerX = sumX / numPoints;
                double centerY = sumY / numPoints;
                double radius = Math.sqrt((sumX2 + sumY2) / numPoints - centerX * centerX - centerY * centerY);

                this.centerX = centerX;
                this.centerY = centerY;
                this.radius = radius;
            }
        }
        public class SpiralRegressionModel extends RegressionModel {
            private double centerX;
            private double centerY;
            private double startRadius;
            private double angularSpeed;
            private double radialSpeed;

            @Override
            public void train(List<Coordinate> coordinates) {
                this.coordinates = coordinates;
                calculateParameters();
            }

            @Override
            public Coordinate predict(int n) {
                if (centerX == 0 && centerY == 0 && startRadius == 0 && angularSpeed == 0 && radialSpeed == 0) {
                    throw new IllegalStateException("Model has not been trained yet.");
                }

                double theta = Math.toRadians(angularSpeed * n);
                double radius = startRadius + radialSpeed * n;
                double x = centerX + radius * Math.cos(theta);
                double y = centerY + radius * Math.sin(theta);

                return new Coordinate(x, y);
            }

            private void calculateParameters() {
                int numPoints = coordinates.size();

                double sumX = 0;
                double sumY = 0;
                double sumX2 = 0;
                double sumY2 = 0;
                double sumXY = 0;
                double sumR = 0;
                double sumR2 = 0;
                double sumRT = 0;

                for (Coordinate coordinate : coordinates) {
                    double x = coordinate.getX();
                    double y = coordinate.getY();
                    double r = Math.sqrt(x * x + y * y);
                    double theta = Math.atan2(y, x);

                    sumX += x;
                    sumY += y;
                    sumX2 += x * x;
                    sumY2 += y * y;
                    sumXY += x * y;
                    sumR += r;
                    sumR2 += r * r;
                    sumRT += r * theta;
                }

                double centerX = sumX / numPoints;
                double centerY = sumY / numPoints;
                double startRadius = sumR / numPoints;
                double angularSpeed = sumRT / sumR2;
                double radialSpeed = Math.sqrt((sumX2 + sumY2) / numPoints - centerX * centerX - centerY * centerY);

                this.centerX = centerX;
                this.centerY = centerY;
                this.startRadius = startRadius;
                this.angularSpeed = angularSpeed;
                this.radialSpeed = radialSpeed;
            }
        }
        public class StillRegressionModel extends RegressionModel {
            private double averageX;
            private double averageY;

            @Override
            public void train(List<Coordinate> coordinates) {
                this.coordinates = coordinates;
                calculateParameters();
            }

            @Override
            public Coordinate predict(int n) {
                if (averageX == 0 && averageY == 0) {
                    throw new IllegalStateException("Model has not been trained yet.");
                }

                return new Coordinate(averageX, averageY);
            }

            private void calculateParameters() {
                int numPoints = coordinates.size();

                double sumX = 0;
                double sumY = 0;

                for (Coordinate coordinate : coordinates) {
                    double x = coordinate.getX();
                    double y = coordinate.getY();

                    sumX += x;
                    sumY += y;
                }

                double averageX = sumX / numPoints;
                double averageY = sumY / numPoints;

                this.averageX = averageX;
                this.averageY = averageY;
            }
        }
    }
}
