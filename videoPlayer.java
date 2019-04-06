package assignmentone;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import static java.lang.Thread.sleep;
import java.util.Arrays;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;

/**
 *
 * @author Arvind Bright
 */
public class videoPlayer {

    JFrame frame;
    JLabel lbIm1;
    BufferedImage imgOne;
    int width = 960;
    int height = 540;
    int FPS;
    double heightScale;
    double widthScale;
    boolean antiAliasing;
    int extraCredit;
    int scaledHeight;
    int scaledWidth;

    private byte[] getRGBfromBytes(byte[] bytes, int x, int y) {
        int r = 0, b = 0, g = 0;
        if (!antiAliasing || x == 0 || x == width - 1 || y == 0 || y == height - 1) {
            return getRGBfromBytes(bytes, x, y, true);
        }

        // anti-aliasing filter
        byte a[];
        for (int i = -1; i <= 1; i++) {
            for (int j = -1; j <= 1; j++) {
                a = getRGBfromBytes(bytes, x + i, y + j, true);
                r += a[0] & 0xff;
                g += a[1] & 0xff;
                b += a[2] & 0xff;
            }
        }
        byte ans[] = {(byte) (r / 9), (byte) (g / 9), (byte) (b / 9)};
        return ans;
    }

    private byte[] getRGBfromBytes(byte[] bytes, int x, int y, boolean recurr) {
        byte r, g, b;
        int offset = width * y + x;
        r = bytes[offset];
        g = bytes[offset + height * width];
        b = bytes[offset + height * width * 2];
        byte a[] = {(byte) r, (byte) g, (byte) b};
        return a;
    }

    /**
     * Read Image RGB Reads the image of given width and height at the given
     * imgPath into the provided BufferedImage.
     */
    private void readImageRGB(String imgPath, BufferedImage img, int frameCount) {

        //to map pixels from file to pixels array with linear scale
        try {
            int frameLength = width * height * 3;

            File file = new File(imgPath);
            RandomAccessFile raf = new RandomAccessFile(file, "r");
            long len = frameLength;
            raf.seek(frameCount * frameLength);
            byte[] bytes = new byte[(int) len];
            raf.read(bytes);
            int Realx, Realy;

            if (extraCredit == 0) {
                for (int y = 0; y < scaledHeight; y++) {
                    for (int x = 0; x < scaledWidth; x++) {
                        byte a[], r, g, b;
                        Realx = (int) Math.floor(((double) width / (double) scaledWidth) * x);
                        Realy = (int) Math.floor(((double) height / (double) scaledHeight) * y);
                        a = getRGBfromBytes(bytes, Realx, Realy);
                        r = a[0];
                        g = a[1];
                        b = a[2];
                        int pix = 0xff000000 | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
                        img.setRGB(x, y, pix);
                    }

                }
            }

            if (extraCredit == 1) {

                double variableSlopeX[] = new double[scaledWidth / 2];
                double variableSlopeY[] = new double[scaledHeight / 2];

                double invHeightScale = 1.0 / heightScale;
                double invWidthScale = 1.0 / widthScale;

                for (int i = 0; i < scaledHeight / 2; i++) {
                    variableSlopeY[i] = invHeightScale + ((1.0 - invHeightScale) / (scaledHeight / 2)) * i;
                }
                for (int i = 0; i < scaledWidth / 2; i++) {
                    variableSlopeX[i] = invWidthScale + ((1.0 - invWidthScale) / (scaledWidth / 2)) * i;
                }

                double realX = 0, realY = 0;

                for (int y = 0; y < scaledHeight; y++) {
                    realX = 0;
                    for (int x = 0; x < scaledWidth; x++) {
                        byte a[], r, g, b;
                        a = getRGBfromBytes(bytes, (int) realX, (int) realY);
                        r = a[0];
                        g = a[1];
                        b = a[2];
                        int pix = 0xff000000 | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
                        img.setRGB(x, y, pix);
                        try {
                            realX += variableSlopeX[x];
                        } catch (Exception e) {
                            realX += variableSlopeX[(scaledWidth) - x - 1];
                        }

                    }
                    try {
                        realY += variableSlopeY[y];
                    } catch (Exception e) {
                        realY += variableSlopeY[(scaledHeight) - y - 1];
                    }

                }

            }

            if (extraCredit == 2) {
                BufferedImage img2 = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
                for (int y = 0; y < height; y++) {

                    for (int x = 0; x < width; x++) {
                        byte a[], r, g, b;
                        a = getRGBfromBytes(bytes, x, y);
                        r = a[0];
                        g = a[1];
                        b = a[2];
                        int pix = 0xff000000 | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
                        img2.setRGB(x, y, pix);
                    }
                }
                int widthDiff = width - scaledWidth;
                int heightDiff = height - scaledHeight;
                SeamCarver sm = new SeamCarver(img2);

                for (int i = 0; i < widthDiff; i++) {
                    sm.removeVerticalSeam(sm.findVerticalSeam());
                }
                for (int i = 0; i < heightDiff; i++) {
                    sm.removeHorizontalSeam(sm.findHorizontalSeam());

                }

                for (int y = 0; y < heightScale; y++) {
                    for (int x = 0; x < widthScale; x++) {
                        img.setRGB(x, y, img.getRGB(x, y));
                    }
                }

            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void showIms(String[] args) {

        // Read a parameter from command line
        String inputFile = args[0];
        widthScale = Double.parseDouble(args[1]);
        heightScale = Double.parseDouble(args[2]);
        FPS = Integer.parseInt(args[3]);
        antiAliasing = (Integer.parseInt(args[4]) == 1);
        extraCredit = (Integer.parseInt(args[5]));

        scaledHeight = (int) Math.round(height * heightScale);
        scaledWidth = (int) Math.round(width * widthScale);
        
        if (extraCredit == 2) {
            extraCreditTwoHelper(inputFile);
            return;
        }
        imgOne = new BufferedImage(scaledWidth, scaledHeight, BufferedImage.TYPE_INT_RGB);
        frame = new JFrame();
        GridBagLayout gLayout = new GridBagLayout();
        frame.getContentPane().setLayout(gLayout);
        GridBagConstraints c = new GridBagConstraints();
        c.fill = GridBagConstraints.HORIZONTAL;
        c.anchor = GridBagConstraints.CENTER;
        c.weightx = 0.5;
        c.gridx = 0;
        c.gridy = 0;

        c.fill = GridBagConstraints.HORIZONTAL;
        c.gridx = 0;
        c.gridy = 1;
        frame.setVisible(true);

        File f = new File(inputFile);
        int numberOfFrames = (int) (f.length() / (height * width * 3));

        for (int i = 0; i < numberOfFrames; i++) {
            long startTime = System.currentTimeMillis();

            readImageRGB(inputFile, imgOne, i);

            // Use label to display the image
            lbIm1 = new JLabel(new ImageIcon(imgOne));

            frame.getContentPane().add(lbIm1, c);

            frame.pack();
            long endTime = System.currentTimeMillis();
            try {
                sleep(((1000 / FPS) - (endTime - startTime)));
            } catch (Exception e) {
                i += Math.round((endTime - startTime) / (1000 / FPS)) - 1;
            }
        }
    }

    private void extraCreditTwoHelper(String inputFile) {
        File f = new File(inputFile);
        int numberOfFrames = (int) (f.length() / (height * width * 3));

        BufferedImage[] img = new BufferedImage[numberOfFrames];

        try {
            frame = new JFrame();
            GridBagLayout gLayout = new GridBagLayout();
            frame.getContentPane().setLayout(gLayout);
            GridBagConstraints c = new GridBagConstraints();
            c.fill = GridBagConstraints.HORIZONTAL;
            c.anchor = GridBagConstraints.CENTER;
            c.weightx = 0.5;
            c.gridx = 0;
            c.gridy = 0;

            c.fill = GridBagConstraints.HORIZONTAL;
            c.gridx = 0;
            c.gridy = 1;

            int frameLength = width * height * 3;

            File file = new File(inputFile);
            RandomAccessFile raf = new RandomAccessFile(file, "r");
            long len = frameLength;
            byte[] bytes = new byte[(int) len];
            System.out.println("please wait until we process all the frames ...");
            System.out.println("Total number of frames: "+numberOfFrames);
            for (int i = 0; i < numberOfFrames; i++) {
                raf.read(bytes);
                img[i] = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);

                for (int y = 0; y < height; y++) {
                    for (int x = 0; x < width; x++) {
                        byte a[], r, g, b;
                        a = getRGBfromBytes(bytes, x, y);
                        r = a[0];
                        g = a[1];
                        b = a[2];
                        int pix = 0xff000000 | ((r & 0xff) << 16) | ((g & 0xff) << 8) | (b & 0xff);
                        img[i].setRGB(x, y, pix);
                    }
                }
                SeamCarver sm = new SeamCarver(img[i]);
                for (int yy = 0; yy < height - scaledHeight; yy++) {
                    sm.removeHorizontalSeam(sm.findHorizontalSeam());

                }
                for (int yy = 0; yy < width - scaledWidth; yy++) {
                    sm.removeVerticalSeam(sm.findVerticalSeam());
                }

                img[i] = sm.picture();

                System.out.println("Read Frames: " + (i+1));
            }
            
            lbIm1 = new JLabel("Test");
            frame.getContentPane().add(lbIm1);

            for (int i = 0; i < numberOfFrames; i++) {
                long startTime = System.currentTimeMillis();

                frame.setVisible(true);
                frame.getContentPane().remove(lbIm1);
                lbIm1 = new JLabel(new ImageIcon(img[i]));

                frame.getContentPane().add(lbIm1, c);

                frame.pack();
                long endTime = System.currentTimeMillis();
                try {
                    sleep(((1000 / FPS) - (endTime - startTime)));
                } catch (Exception e) {
                }
            }

        } catch (Exception e) {
            System.err.println(e);
        }

    }

    public static void main(String[] args) {
        AssignmentOne ren = new AssignmentOne();
        ren.showIms(args);
        System.exit(0);
    }

}

class SeamCarver {

    // The representation of the given image
    private int[][] color;

    // The energy of each pixel in the image
    private double[][] energy;

    // Arrays and sinks for finding the shortest path through the image energy
    private double[][] distTo;
    private double distToSink;
    private int[][] edgeTo;
    private int edgeToSink;

    // The current width and height
    private int w;
    private int h;

    // False if finding or removing a vertical seam,
    // true if finding or removing a horizontal seam.
    private boolean transposed;

    /**
     * Create a seam carver object based on the given picture.
     *
     * @param picture the given picture
     * @throws NullPointerException if the given picture is {@code null}.
     */
    public SeamCarver(BufferedImage picture) {
        if (picture == null) {
            throw new java.lang.NullPointerException();
        }

        // Initialize the dimensions of the picture
        w = picture.getWidth();
        h = picture.getHeight();

        // Store the picture's color information in an int array,
        // using the RGB coding described at:
        // http://docs.oracle.com/javase/8/docs/api/java/awt/Color.html#getRGB()
        color = new int[h][w];

        // Set the dimensions of the energy array
        energy = new double[h][w];

        // Store color information
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                color[i][j] = picture.getRGB(j, i);
            }
        }

        // Pre-calculate the energy array
        for (int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                energy[i][j] = calcEnergy(j, i);
            }
        }
    }

    /**
     * Current picture.
     *
     * @return the current picture.
     */
    public BufferedImage picture() {

        // Create and return a new pic with the stored color information
        BufferedImage pic = new BufferedImage(width(), height(), BufferedImage.TYPE_INT_RGB);
        for (int i = 0; i < height(); i++) {
            for (int j = 0; j < width(); j++) {
                pic.setRGB(j, i, new Color(color[i][j]).getRGB());
            }
        }
        return pic;
    }

    /**
     * Width of current picture.
     *
     * @return the width of the current picture.
     */
    public int width() {
        return w;
    }

    /**
     * Height of current picture.
     *
     * @return the height of the current picture.
     */
    public int height() {
        return h;
    }

    /**
     * Energy of pixel at column x and row y.
     *
     * Note that (0, 0) is the pixel at the top-left corner of the image.
     *
     * The dual-gradient energy function is used to compute the energy of a
     * pixel.
     *
     * @param x
     * @param y
     * @return the energy of the pixel at column <em>x</em> and row <em>y</em>.
     * @throws IndexOutOfBoundsException if <em>x</em> is greater than or equal
     * to the image width, if <em>y</em> is greater than or equal to the image
     * height, or if <em>x</em> or <em>y</em> are negative.
     */
    public double energy(int x, int y) {
        if (x >= width() || y >= height() || x < 0 || y < 0) {
            throw new java.lang.IndexOutOfBoundsException();
        }

        return energy[y][x];
    }

    /**
     * Helper method to calculate the energy of pixel at column x and row y.
     *
     * Note that (0,0) is the pixel at the top-left corner of the image.
     *
     * The dual-gradient energy function is used to compute the energy of a
     * pixel.
     *
     * @param x
     * @param y
     * @return the energy of the pixel at column <em>x</em> and row <em>y</em>.
     * @throws IndexOutOfBoundsException if <em>x</em> is greater than or equal
     * to the image width, if <em>y</em> is greater than or equal to the image
     * height, or if <em>x</em> or <em>y</em> are negative.
     */
    private double calcEnergy(int x, int y) {
        if (x >= width() || y >= height() || x < 0 || y < 0) {
            throw new java.lang.IndexOutOfBoundsException();
        }

        // Return 1000.0 for border pixels
        if (x == 0 || y == 0 || x == width() - 1 || y == height() - 1) {
            return (double) 1000;
        }

        // Store pixel values in Color objects.
        Color up = new Color(color[y - 1][x]);
        Color down = new Color(color[y + 1][x]);
        Color left = new Color(color[y][x - 1]);
        Color right = new Color(color[y][x + 1]);

        return Math.sqrt(gradient(up, down) + gradient(left, right));
    }

    /**
     * Returns the gradient computed from the two Colors <em>a</em> and
     * <em>b</em>.
     *
     * @param a the first Color
     * @param b the second Color
     * @return the gradient of <em>a</em> and <em>b</em>.
     */
    private double gradient(Color a, Color b) {
        return Math.pow(a.getRed() - b.getRed(), 2)
                + Math.pow(a.getBlue() - b.getBlue(), 2)
                + Math.pow(a.getGreen() - b.getGreen(), 2);
    }

    /**
     * Sequence of indices for horizontal seam.
     *
     * This method conducts a shortest-path search as if the energy matrix were
     * an edge-weighted directed acyclic graph.
     *
     * The source vertex is an implicit vertex sitting to the left of the image,
     * to which all of the left-column pixels are adjacent.
     *
     * The sink vertex is an explicit vertex sitting to the right of the image,
     * which is (the only vertex) adjacent to all of the right-column pixels.
     *
     * Each pixel can visit only the pixel to its immediate right, the pixel to
     * its right and above it (if possible), and the pixel to its right and
     * below it (if possible).
     *
     * @return the sequence of indices for the horizontal seam.
     */
    public int[] findHorizontalSeam() {
        transposed = true;

        // Reset our distTo and edgeTo values for a new search
        distToSink = Double.POSITIVE_INFINITY;
        edgeToSink = Integer.MAX_VALUE;
        distTo = new double[h][w];
        edgeTo = new int[h][w];
        for (double[] r : distTo) {
            Arrays.fill(r, Double.POSITIVE_INFINITY);
        }
        for (int[] r : edgeTo) {
            Arrays.fill(r, Integer.MAX_VALUE);
        }

        // Relax the entire left column, since this is our starting column
        for (int i = 0; i < height(); i++) {
            distTo[i][0] = (double) 1000;
            edgeTo[i][0] = -1;
        }

        // Visit all pixels from the left side, diagonally to the right,
        // in keeping with topological order.
        // The topological order is the reverse of the DFS post-order,
        // which visits the top-most adjacent pixel first, before it visits
        // the pixels below.
        for (int depth = height() - 1; depth > 0; depth--) {
            for (int out = 0;
                    out < width() && depth + out < height();
                    out++) {
                visit(depth + out, out);
            }
        }

        // Visit all pixels from the top, diagonally to the right,
        // in keeping with the topological order described above.
        for (int top = 0; top < width(); top++) {
            for (int depth = 0;
                    depth + top < width() && depth < height();
                    depth++) {
                visit(depth, depth + top);
            }
        }

        // Populate seam[] with the shortest path
        int[] seam = new int[width()];
        seam[width() - 1] = edgeToSink;

        for (int j = width() - 1; j > 0; j--) {
            seam[j - 1] = edgeTo[seam[j]][j];
        }

        // null out our shortest-path arrays for garbage collection
        distTo = null;
        edgeTo = null;

        return seam;

    }

    /**
     * Sequence of indices for vertical seam.
     *
     * This method conducts a shortest-path search as if the energy matrix were
     * an edge-weighted directed acyclic graph.
     *
     * The source vertex is an implicit vertex sitting above the image, to which
     * all of the top-row pixels are adjacent.
     *
     * The sink vertex is an explicit vertex sitting below the image, which is
     * (the only vertex) adjacent to all of the bottom-row pixels.
     *
     * Each pixel can visit only the pixel directly below it, the pixel below it
     * and to its left (if possible), and the pixel below it and to its right
     * (if possible).
     *
     * @return the sequence of indices for the vertical seam.
     */
    public int[] findVerticalSeam() {
        transposed = false;

        // Reset our distTo and edgeTo values for a new search
        distToSink = Double.POSITIVE_INFINITY;
        edgeToSink = Integer.MAX_VALUE;
        distTo = new double[h][w];
        edgeTo = new int[h][w];
        for (double[] r : distTo) {
            Arrays.fill(r, Double.POSITIVE_INFINITY);
        }
        for (int[] r : edgeTo) {
            Arrays.fill(r, Integer.MAX_VALUE);
        }

        // Relax the entire top row, since this is our starting row
        Arrays.fill(distTo[0], (double) 1000);
        Arrays.fill(edgeTo[0], -1);

        // Visit all pixels from the top, diagonally to the right,
        // in keeping with topological order.
        // The topological order is the reverse of the DFS post-order,
        // which visits the left-most adjacent pixel first, before it visits
        // pixels to the right.
        for (int top = width() - 1; top >= 0; top--) {
            for (int depth = 0;
                    depth + top < width() && depth < height();
                    depth++) {
                visit(depth, depth + top);
            }
        }
        // Visit all pixels from the left side, diagonally to the right,
        // in keeping with the topological order described above.
        for (int depth = 1; depth < height(); depth++) {
            for (int out = 0;
                    out < width() && depth + out < height();
                    out++) {
                visit(depth + out, out);
            }
        }

        // Populate seam[] with the shortest path
        int[] seam = new int[height()];
        seam[height() - 1] = edgeToSink;

        for (int i = height() - 1; i > 0; i--) {
            seam[i - 1] = edgeTo[i][seam[i]];
        }

        // null out our shortest-path arrays for garbage collection
        distTo = null;
        edgeTo = null;

        return seam;
    }

    /**
     * Given a pixel's coordinates, relax the pixels adjacent to that pixel.
     *
     * @param i the vertical index of the pixel
     * @param j the horizontal index of the pixel
     */
    private void visit(int i, int j) {
        if (transposed) {
            // Only relax the sink
            if (j == width() - 1) {
                relax(i, j);
            } // Bottom edge; relax to the right and above
            else if (i == height() - 1) {
                relax(i, j, i, j + 1);
                relax(i, j, i - 1, j + 1);
            } // Top edge; relax to the right and below
            else if (i == 0) {
                relax(i, j, i, j + 1);
                relax(i, j, i + 1, j + 1);
            } // Middle pixel; relax right, below, and above
            else {
                relax(i, j, i - 1, j + 1);
                relax(i, j, i, j + 1);
                relax(i, j, i + 1, j + 1);
            }
        } else {
            // Only relax the sink
            if (i == height() - 1) {
                relax(i, j);
            } // Right edge; relax below and to the left
            else if (j == width() - 1) {
                relax(i, j, i + 1, j - 1);
                relax(i, j, i + 1, j);
            } // Left edge; relax below and to the right
            else if (j == 0) {
                relax(i, j, i + 1, j);
                relax(i, j, i + 1, j + 1);
            } // Middle pixel; relax left, below, and right
            else {
                relax(i, j, i + 1, j - 1);
                relax(i, j, i + 1, j);
                relax(i, j, i + 1, j + 1);
            }
        }
    }

    /**
     * Given an index, relax the sink vertex from that index.
     *
     * This method should only be called on the "last" vertices in the image.
     *
     * @param i the vertical index of the pixel
     * @param j the horizontal index of the pixel
     */
    private void relax(int i, int j) {
        if (validIndex(i, j)) {
            if (distToSink > distTo[i][j]) {
                distToSink = distTo[i][j];
                if (transposed) {
                    edgeToSink = i;
                } else {
                    edgeToSink = j;
                }
            }
        }
    }

    /**
     * Given the coordinates of pixel 1 and pixel 2, relax pixel 2 from pixel 1.
     *
     * This method should not be called on the "last" pixels in the image.
     *
     * @param i1 the vertical index of pixel 1
     * @param j1 the horizontal index of pixel 1
     * @param i2 the vertical index of pixel 2
     * @param j2 the horizontal index of pixel 2
     */
    private void relax(int i1, int j1, int i2, int j2) {
        if (validIndex(i1, j1) && validIndex(i2, j2)) {
            if (distTo[i2][j2] > distTo[i1][j1] + energy[i2][j2]) {
                distTo[i2][j2] = distTo[i1][j1] + energy[i2][j2];
                if (transposed) {
                    edgeTo[i2][j2] = i1;
                } else {
                    edgeTo[i2][j2] = j1;
                }
            }
        }
    }

    /**
     * Is the given pixel coordinate valid?
     *
     * @param i the vertical index of the pixel
     * @param j the horizontal index of the pixel
     * @return {@code true} if the current picture contains the pixel
     * coordinate, {@code false} otherwise.
     */
    private boolean validIndex(int i, int j) {
        return (i >= 0 && i < height() && j >= 0 && j < width());
    }

    /**
     * Remove horizontal seam from current picture.
     *
     * @param seam the given seam.
     * @throws NullPointerException if the given <em>seam</em> is {@code null}.
     * @throws IllegalArgumentException if the given <em>seam</em> does not
     * match the picture width, if an index in the given <em>seam</em>
     * is negative or is taller than the picture, or if two adjacent entries in
     * the given <em>seam</em> differ by more than 1.
     */
    public void removeHorizontalSeam(int[] seam) {

        // Check for bad input
        if (height() <= 1) {
            throw new java.lang.IllegalArgumentException("Picture too short");
        }
        if (seam == null) {
            throw new java.lang.NullPointerException();
        }
        if (seam.length != width()) {
            throw new java.lang.IllegalArgumentException("Invalid seam length");
        }

        int yLast = seam[0];
        for (int y : seam) {
            if (y >= height() || y < 0) {
                throw new java.lang.IllegalArgumentException("Index out of bounds");
            }
            if (Math.abs(y - yLast) > 1) {
                throw new java.lang.IllegalArgumentException("Index not adjacent");
            }
            yLast = y;
        }

        // Create replacement arrays
        int[][] newColor = new int[height() - 1][width()];
        double[][] newEnergy = new double[height() - 1][width()];

        // Populate replacement arrays, skipping pixels in the seam
        for (int j = 0; j < width(); j++) {
            int s = seam[j];
            for (int i = 0; i < s; i++) {
                newColor[i][j] = color[i][j];
                newEnergy[i][j] = energy[i][j];
            }

            for (int i = s + 1; i < height(); i++) {
                newColor[i - 1][j] = color[i][j];
                newEnergy[i - 1][j] = energy[i][j];
            }
        }

        color = newColor;
        energy = newEnergy;
        h--;

        // Recalculate the energy along the seam
        for (int j = 0; j < width(); j++) {
            int s = seam[j];
            // Top edge removed
            if (s == 0) {
                energy[s][j] = calcEnergy(j, s);
            } // Bottom edge removed
            else if (s == height()) {
                energy[s - 1][j] = calcEnergy(j, s - 1);
            } // Middle pixel removed
            else {
                energy[s][j] = calcEnergy(j, s);
                energy[s - 1][j] = calcEnergy(j, s - 1);
            }
        }
    }

    /**
     * Remove vertical seam from current picture.
     *
     * @param seam the given seam.
     * @throws NullPointerException if the given <em>seam</em> is {@code null}.
     * @throws IllegalArgumentException if the given <em>seam</em> does not
     * match the picture height, if an index in the given <em>seam</em>
     * is negative or is wider than the picture, or if two adjacent entries in
     * the given <em>seam</em> differ by more than 1.
     */
    public void removeVerticalSeam(int[] seam) {

        // Check for bad input
        if (width() <= 1) {
            throw new java.lang.IllegalArgumentException("Picture too narrow");
        }
        if (seam == null) {
            throw new java.lang.NullPointerException();
        }
        if (seam.length != height()) {
            throw new java.lang.IllegalArgumentException("Invalid seam length");
        }

        int xLast = seam[0];
        for (int x : seam) {
            if (x >= width() || x < 0) {
                throw new java.lang.IllegalArgumentException("Index out of bounds");
            }
            if (Math.abs(x - xLast) > 1) {
                throw new java.lang.IllegalArgumentException("Index not adjacent");
            }
            xLast = x;
        }

        // Create replacement arrays
        int[][] newColor = new int[height()][width() - 1];
        double[][] newEnergy = new double[height()][width() - 1];

        // Populate replacement arrays, skipping pixels in the seam
        for (int i = 0; i < height(); i++) {
            int s = seam[i];

            for (int j = 0; j < s; j++) {
                newColor[i][j] = color[i][j];
                newEnergy[i][j] = energy[i][j];
            }

            for (int j = s + 1; j < width(); j++) {
                newColor[i][j - 1] = color[i][j];
                newEnergy[i][j - 1] = energy[i][j];
            }
        }

        color = newColor;
        energy = newEnergy;
        w--;

        // Recalculate the energy along the seam
        for (int i = 0; i < height(); i++) {
            int s = seam[i];

            // Left edge removed
            if (s == 0) {
                energy[i][s] = calcEnergy(s, i);
            } // Right edge removed
            else if (s == width()) {
                energy[i][s - 1] = calcEnergy(s - 1, i);
            } // Middle pixel removed
            else {
                energy[i][s] = calcEnergy(s, i);
                energy[i][s - 1] = calcEnergy(s - 1, i);
            }
        }
    }

}
