#include <stdio.h>
#include <stdlib.h>
//Change these as per the requirement of the project datasets
#define rows 7129
#define cols 34
int main()
{
    // Open the text file for reading
    FILE *fp_in = fopen("training.txt", "r");
    if (fp_in == NULL) {
        printf("Error: Unable to open the input file\n");
        return 1;
    }

    // Open the binary file for writing
    FILE *fp_out = fopen("output.bin", "wb");
    if (fp_out == NULL) {
        printf("Error: Unable to open the output file\n");
        return 1;
    }

    // Read and write the data
    // The data of the genes we read in is of type double 
    double data;
 
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if (fscanf(fp_in, "%lf", &data) == 1) {
                fwrite(&data, sizeof(double), 1, fp_out);
            } else {
                printf("Error: Unable to read the input file\n");
                fclose(fp_in);
                fclose(fp_out);
              exit(-1);
            }
        }
    }

    // Close the files
    fclose(fp_in);
    fclose(fp_out);

    printf("Binary file is created!\n");

    return 0;
}
