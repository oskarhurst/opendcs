import java.io.File;  // Import the File class
import java.io.FileWriter;   // Import the FileWriter class
import java.io.IOException;  // Import the IOException class to handle errors
import decodes.tsdb.comprungui.CompExecTest;

public class LocalAlgoTest {
    public static void main(String[] args) {
//        try {
//            String pathname = "filename.txt";
//            File myObj = new File(pathname);
//            myObj.createNewFile();
//            FileWriter myWriter = new FileWriter(pathname);
//            myWriter.write("Files in Java might be tricky, but it is fun enough!");
//            myWriter.close();
//            System.out.println("Successfully wrote to the file.");
//            if (myObj.delete()) {
//                System.out.println("Deleted the file: " + myObj.getName());
//            } else {
//                System.out.println("Failed to delete the file.");
//            }
//        } catch (IOException e) {
//            System.out.println("An error occurred.");
//            e.printStackTrace();
//        }

        CompExecTest.main(args);
    }
}
