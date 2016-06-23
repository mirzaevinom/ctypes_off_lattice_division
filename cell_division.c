#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>



//Used "Valgrind" to find memory leaks in the code

// This function checks for empty spots around the cell undergoing mitosis

int empty_spot_check( int a_size, //location matrix size
              int b_size, //size of spheric shift matrix
              int cnum,   //cell number in mitosis
              double temp_rad, //trial radius
              double old_rad, //old radius of the cell going through mitosis
              int *arr_shuffle, //stores the index of empty spotss
              double loc_mat[][5], //location matrix should be given as an input
              double **sphr_shift //spherical shift matrix
             )
{
  int i , j, esptnum=0, dir_list;
  double espots[3];

  for (i = 0; i < b_size; i++)
    {
		dir_list=1;

		//Initialization of empty spots around the cell

		espots[0] = loc_mat[cnum][0] + (3 * temp_rad - old_rad) * sphr_shift[i][0];
        espots[1] = loc_mat[cnum][1] + (3 * temp_rad - old_rad) * sphr_shift[i][1];
        espots[2] = loc_mat[cnum][2] + (3 * temp_rad - old_rad) * sphr_shift[i][2];


        for (j=0; j<a_size; j++)
        {
            dir_list = dir_list * ( sqrt ( (loc_mat[j][0]- espots[0])*(loc_mat[j][0]- espots[0]) + (loc_mat[j][1] - espots[1])*(loc_mat[j][1] - espots[1]) + (loc_mat[j][2] - espots[2])*(loc_mat[j][2] - espots[2]) ) >= ( temp_rad  + loc_mat[j][3] ) );
        }

       //If the i-th spot around the cell is empty, the i'th index will be stored in the arr_shuffle index
       if ( dir_list == 1 )
       {
		   arr_shuffle[esptnum]=i;
		   esptnum++;
	   }

    }
    //Returns number of empty spots around the cell
    return esptnum;

}

//Given an array with n elements, the function returns shuffled array.
//Randomization depends on the time of the day.
//The function is modified version of the one on stackexchange.
// See actual code here, http://stackoverflow.com/questions/6127503/shuffle-array-in-c

void random_shuffle(int *array, size_t n)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    int usec = tv.tv_usec;
    srand48(usec);


    if (n > 1) {
        size_t i;
        for (i = n - 1; i > 0; i--) {
            size_t j = (unsigned int) (drand48()*(i+1));
            int t = array[j];
            array[j] = array[i];
            array[i] = t;
        }
    }
}


//This function returns total volume of the cells at each step

double vol_sum(double a[][5], int rows)
{
   int i;
   double sum=0;
   for (i=0; i<rows; i++)
   {

	 sum = sum + pow( a[i][3] , 3);

   }
   return(sum);
}



int cell_division( int sizeofrad, // length of allowed radius vector
                   int sizeofloc, // maximum length that location matrix can grow
                   int numcells,  // number of cells in the beginning
                   int numsteps, //number of steps the simulation run
                   int num_phi,  //number of discretization points for angle phi of spherical coordinate
                   int num_theta, //number of discretization points for angle theta of spherical coordinate
                   double loc_mat[][5], //pointer for location matrix passed by python code
                   double allowed_rad[], //list of allowed radius for daughter cell
                   double volume[] //stores the total volume of the cells at each step
                   )
{
    int steps, cnum , numrad , esp=0 , trvs, numiter, i, espt_number;
    int sizeofshift = num_phi*num_theta;

    double temp_rad, old_rad;

    int * arr_shuffle;
    arr_shuffle = (int *) malloc(sizeofshift * sizeof(int));

    for (i=0; i<sizeofshift; i++)
        arr_shuffle[i]=i;

	int *cell_shuffle;
	cell_shuffle = (int * ) malloc(sizeofloc * sizeof(int) );


    double *phi , *theta ;
    phi = (double *) malloc(num_phi*num_theta * sizeof(double));
    theta = (double *) malloc(num_phi*num_theta * sizeof(double));


    double **sphr_shift = (double **) malloc(num_phi*num_theta * sizeof(double *));

	for(i = 0; i < num_phi*num_theta; i++)
		sphr_shift[i] = (double * ) malloc(3 * sizeof(double));


	//initializes spherical shift matrix

    for (i=0; i<num_phi*num_theta; i++)
    {

        phi[i] = M_PI / num_phi * (i % num_phi);

        theta[i] = 2*M_PI / num_theta  * (i % num_theta);
        sphr_shift[i][0] = cos(theta[i]) * sin(phi[i]);
        sphr_shift[i][1] = sin(theta[i]) * sin(phi[i]);
        sphr_shift[i][2] = cos(phi[i]);
    }






    for (steps=0; steps<numsteps; steps++)
    {
        //Sets the number of iteration for each step to the number of cells from the previous step
        numiter=numcells;

        //An array used to randomly access existing cells
        for (i=0; i<numiter; i++)
			cell_shuffle[i]=i;

        random_shuffle(cell_shuffle, numiter);

        for (i=0; i<numiter; i++)
        {
			cnum = cell_shuffle[i];

            if (loc_mat[cnum][4]==1 && numcells< (sizeofloc-1) )
            {
				//Empty spots around the cells are checked for each allowed radius
                for (numrad=0; numrad<sizeofrad; numrad++)
                {

                    temp_rad = allowed_rad[numrad];
                    old_rad  = loc_mat[cnum][3];

                    //To avoid empty_spot_check radius of the cell in mitosis is set to zero
                    loc_mat[cnum][3] = 0;


                    // This function checks for empty spots around the cell undergoing mitosis


                    int i , j, esptnum=0, dir_list;
					double espots[3];

					for (i = 0; i < sizeofshift; i++)
						{
							dir_list=1;

							//Initialization of empty spots around the cell

							espots[0] = loc_mat[cnum][0] + (3 * temp_rad - old_rad) * sphr_shift[i][0];
							espots[1] = loc_mat[cnum][1] + (3 * temp_rad - old_rad) * sphr_shift[i][1];
							espots[2] = loc_mat[cnum][2] + (3 * temp_rad - old_rad) * sphr_shift[i][2];


							for (j=0; j<numcells; j++)
								{
									dir_list = dir_list * ( sqrt ( (loc_mat[j][0]- espots[0])*(loc_mat[j][0]- espots[0]) + (loc_mat[j][1] - espots[1])*(loc_mat[j][1] - espots[1]) + (loc_mat[j][2] - espots[2])*(loc_mat[j][2] - espots[2]) ) >= ( temp_rad  + loc_mat[j][3] ) );
								}

							//If the i-th spot around the cell is empty, the i'th index will be stored in the arr_shuffle index
							if ( dir_list == 1 )
								{
									arr_shuffle[esptnum]=i;
									esptnum++;
								}

						}



					//if the number of empty spots around the cell is not equal to zero, the code puts new cell in one of the empty spots

                    if ( espt_number != 0  )
                       {


                           //random_shuffle( arr_shuffle, espt_number );

                           esp= rand() % espt_number;

                           loc_mat[ numcells ][0] = loc_mat[cnum][0] + (3 * temp_rad - old_rad) * sphr_shift[ arr_shuffle[esp] ][0];
                           loc_mat[ numcells ][1] = loc_mat[cnum][1] + (3 * temp_rad - old_rad) * sphr_shift[ arr_shuffle[esp] ][1];
                           loc_mat[ numcells ][2] = loc_mat[cnum][2] + (3 * temp_rad - old_rad) * sphr_shift[ arr_shuffle[esp] ][2];
                           loc_mat[ numcells ][3] = temp_rad;
                           loc_mat[ numcells ][4] = 1;
                           numcells++;

                           loc_mat[ cnum ][3] = temp_rad;
                           loc_mat[ cnum ][0] = loc_mat[ cnum ][0] - ( old_rad - temp_rad) * sphr_shift[ arr_shuffle[esp] ][0];
                           loc_mat[ cnum ][1] = loc_mat[ cnum ][1] - ( old_rad - temp_rad) * sphr_shift[ arr_shuffle[esp] ][1];
                           loc_mat[ cnum ][2] = loc_mat[ cnum ][2] - ( old_rad - temp_rad) * sphr_shift[ arr_shuffle[esp] ][2];


                           break;
                       }
                    else
                       {
                          //if there is no empty spot for this allowed radius, radius of the cell in mitosis is set back to original radius
                           loc_mat[cnum][3] = old_rad;

                       }


                    if ((temp_rad==allowed_rad[sizeofrad-1]) && espt_number == 0 )
                        {
							//if there is no empty spot for each allowed radius, the cell enters quiscent state
                            loc_mat[cnum][4] = 0;
                            break;
                        }

                }


            }
        }

        //Total volume of the cells is calculated
        volume[steps] = 4 * M_PI / 3 * vol_sum(loc_mat , numcells);
    }

    free(arr_shuffle);
    free(cell_shuffle);

    for(i = 0; i < num_phi*num_theta; i++)
		free(sphr_shift[i]);

	free(sphr_shift);

	free(phi);
	free(theta);

	//Returns number of cells at the end
    return numcells;

}
