#include "hdf5Builder.h"
#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

const H5std_string FILE_NAME("testing.h5");
const H5std_string DATASET_NAME ( " testArray" );

int hdf5Test(void){
    /*
    * Data initialization.
    */
    int i, j;
    int data[5][5]; // buffer for data to write
    for (j = 0; j < 5; j++){
        for (i = 0; i < 5; i++){
            data[j][i] = i + j;
        }
    }

    // Try block to detect exceptions raised by any of the calls inside it
    try
    {
    /*
    * Turn off the auto-printing when failure occurs so that we can
    * handle the errors appropriately
    */
    Exception::dontPrint();
    /*
    * Create a new file using H5F_ACC_TRUNC access,
    * default file creation properties, and default file
    * access properties.
    */
    H5File file( FILE_NAME, H5F_ACC_TRUNC );
    /*
    * Define the size of the array and create the data space for fixed
    * size dataset.
    */
    hsize_t dimsf[2]; // dataset dimensions
    dimsf[0] = 5;
    dimsf[1] = 5;
    DataSpace dataspace( 2, dimsf );
    /*
    * Define datatype for the data in the file.
    * We will store little endian INT numbers.
    */
    IntType datatype( PredType::NATIVE_INT );
    datatype.setOrder( H5T_ORDER_LE );
    /*
    * Create a new dataset within the file using defined dataspace and
    * datatype and default dataset creation properties.
    */
    DataSet dataset = file.createDataSet( DATASET_NAME, datatype, dataspace );
    /*
    * Write the data to the dataset using default memory space, file
    * space, and transfer properties.
    */
    dataset.write( data, PredType::NATIVE_INT );
    } // end of try block
    // catch failure caused by the H5File operations
    catch( FileIException error )
    {
    error.printError();
    return -1;
    }
    // catch failure caused by the DataSet operations
    catch( DataSetIException error )
    {
    error.printError();
    return -1;
    }
    // catch failure caused by the DataSpace operations
    catch( DataSpaceIException error )
    {
    error.printError();
    return -1;
    }
    // catch failure caused by the DataSpace operations
    catch( DataTypeIException error )
    {
    error.printError();
    return -1;
    }
    return 0; // successfully terminated
}
