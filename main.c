#include <stdio.h>
#include <stdlib.h>

    void drawArr(float** arr, size_t m, size_t n)
    {
        size_t i, j;
        printf("*       B\\/     ");

        for (j = 2; j < n; j++)
            printf("X%1.0f      ", arr[0][j]);

        printf("\nF>   ");
        for (j = 1; j < n; j++)
            printf("%7.3f ", arr[1][j]);
        printf("\n");

        for (i = 2; i < m; i++) {
                printf("X%1.0f   ", arr[i][0]);

            for(j = 1; j < n; j++)
                printf("%7.3f ", arr[i][j]);

            printf("\n");
        }
        printf("\n");
    }

    float** initMemArr(size_t m, size_t n)
    {
        size_t i;
        float** arr = (float**) malloc(m * sizeof(float*));
        for(i = 0; i < m; i++)
            arr[i] = (float*) malloc(n * sizeof(float));

        return arr;
    }

    float** reInitMemArr(float** arr, size_t newM, size_t newN, size_t oldM, size_t oldN)
    {
        size_t i;

        if (oldM != newM) {
            arr = (float**) realloc(arr, newM * sizeof(float*));

            for(i = oldM; i < newM; i++)
                arr[i] = (float*) malloc(newN * sizeof(float));
        }

        if (oldN != newN)
            for(i = 0; i < oldM; i++)
                arr[i] = (float*) realloc(arr[i], newN * sizeof(float));

        return arr;
    }

    void freeMemArr(float** arr, size_t m)
    {
        size_t i;
        for(i = 0; i < m; i++)
            free(arr[i]);
        free(arr);
    }

    /* in-function dynamic array initialization */
    void initValArr(float** trgArr, size_t m, size_t n, float srcArr[m][n])
    {
        size_t i, j;
        for(i = 0; i < m; i++)
            for(j = 0; j < n; j++)
                trgArr[i][j] = srcArr[i][j];
    }

    void copyValArr(float** srcArr, float** trgArr, size_t m, size_t n)
    {
        size_t i, j;
        for(i = 0; i < m; i++)
            for(j = 0; j < n; j++)
                trgArr[i][j] = srcArr[i][j];
    }

    /* ---------------------------------------- */

    void simplexCalc(float** table, size_t m, size_t n, int toMax, int invIndexColumnSigns)
    {
        float** tmpTable,
            un, tmpUn;
        size_t i, j, kI, kJ;
        int isOpt, firstNoCheck;

        tmpTable = initMemArr(m, n);
        /* we think, is not optimized on start part */
        isOpt = 0;

        printf("Simplex calculation by %s extremum.\n", toMax? "max" : "min");

        if (invIndexColumnSigns) {
            for(j = 2; j < m; j++)
                table[1][j] = -table[1][j];

            if (table[1][j] == -0)
                table[1][j] = 0;
        }

        for(;;) {

            drawArr(table, m, n);

            /* Getting key column and checking table by optimization aspect*/
            if (toMax) {
                un = 1;
                for(j = 2; j < n; j++)
                    if (table[1][j] < 0 && table[1][j] < un) {
                        un = table[1][j];
                        kJ = j;
                    }
                /* is Optimal */
                if (un == 1)
                    isOpt = 1;

            /* is to min */
            } else {

                un = -1;
                for(j = 2; j < n; j++)
                    if (table[1][j] > 0 && table[1][j] > un) {
                        un = table[1][j];
                        kJ = j;
                    }
                if (un == -1)
                    isOpt = 1;
            }

        /* Print case */
        for(i = 2; i < m; i++)
            printf("X%.0f = %3.3f, ", table[i][0], table[i][1]);
        printf("\n");
        for(j = 2; j < n; j++)
            printf("X%.0f = 0, ", table[0][j]);
        printf("F(X) = %3.3f\n", table[1][1]);
        /*if its optimized, then finish from simplex-method */
        if (isOpt) {
            printf("Its optimized case plan\n");
            break;
        } else
            printf("Its not optimized case plan\n");

         /* Getting key row */

         firstNoCheck = 1;

         for(i = 2; i < n; i++) {
            tmpUn = table[i][1] / table[i][kJ];
            if (tmpUn < 0)
                tmpUn = -tmpUn;

            if (firstNoCheck || tmpUn < un) {
                un = tmpUn;
                kI = i;
                firstNoCheck = 0;
            }
         }

        printf("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");
        printf("key i: %d, key j: %d\nkey value %3.6f\n", kI - 1, kJ - 1, table[kI][kJ]);
        printf("-----------------------------\n");

        /* Creating new simplex table */
        /*swap variable labels by key indexes*/

        tmpUn = table[0][kJ];
        table[0][kJ] = table[kI][0];
        table[kI][0] = tmpUn;

        for (i = 1; i < m; i++)
            for(j = 1; j < n; j++) {
                if (i == kI && j == kJ)
                    tmpTable[i][j] = 1 / table[i][j];
                else if (i == kI)
                    tmpTable[i][j] = table[i][j] / table[i][kJ];
                else if (j == kJ)
                    tmpTable[i][j] = -(table[i][j] / table[kI][j]);
                else
                    tmpTable[i][j] = (table[i][j] * table[kI][kJ] - table[kI][j] * table[i][kJ]) / table[kI][kJ];

                if (tmpTable[i][j] == -0)
                    tmpTable[i][j] = 0;
            }

         /* Move table on main table var */
            for (i = 1; i < m; i++)
                for(j = 1; j < n; j++)
                    table[i][j] = tmpTable[i][j];

        }   /* End of iteration */
        freeMemArr(tmpTable, m);
    }

    void gomoriCalc(float** table, size_t m, size_t n, int toMax, size_t numVars)
    {
        size_t i, mfI, j;
        int isCase, invertDS;
        float fract, maxFract;
// 1

        if (toMax)
            invertDS = 1;
        else
            invertDS = 0;

        simplexCalc(table, m, n, toMax, 0);

        for(;;) {
// 2
            printf("\nGomori method part\nChecking validation of solution\n");

            isCase = 1;
            maxFract = 0.0f;

            for(i = 2; i < m; i++) {
                fract = table[i][1] - (float)((int)table[i][1]);

                if (fract) {
                    isCase = 0;
                    if (maxFract < fract) {
                        maxFract = fract;
                        mfI = i;
                    }

                    printf("%3.3f is fractional\n", table[i][1]);

                }   else

                printf("%3.1f is integer\n", table[i][1]);
            }
// 3:1
            if (isCase) {
                printf("Is valid solution\n");
                break;  // out of Gomori loop
            }
// 3:2
            else {
                printf("Is not valid solution\nMaximum fract: %0.3f (of %3.3f), i: %d (X%1.0f)\n", maxFract, table[mfI][1], mfI, table[mfI][0]);

                printf("Creating of new constraint\n");

                i = m;
                m++;
                table = reInitMemArr(table, m, n, i, n);    // i is old rows number

                table[i][0] = ++numVars;    // set new variable name (id)

                for (j = 1; j < n; j++) {
                    fract = table[mfI][j] - (float)( (int)table[mfI][j] );
                    if (fract > 0)
                        table[i][j] = -fract;
                    else if (fract < 0)
                        table[i][j] = -(1.0f + fract);
                    else
                        table[i][j] = 0;
                }

                drawArr(table, m, n);

            }
            simplexCalc(table, m, n, 0, invertDS); // to min and invert if its max by once stage
            if (invertDS)
                invertDS = 0;
        }

    }

int main()
{
    float** table;
    size_t m, n;

    m = 4;
    n = 4;
    table = initMemArr(m, n);

    /*
    initValArr(table, m, n,
        (float[4][4]) {
            {0, 0, 1, 2},
            {0, 0, -3, -4},
            {3, 10, 1, 4},
            {4, 8, 3, 2}
        }
    );
*/

    initValArr(table, m, n,
        (float[4][4]) {
            {0, 0, 1, 2},
            {0, 0, -3, -4},
            {3, 10, 1, 4},
            {4, 8, 3, 2}
        }
    );


    gomoriCalc(table, m, n, 1, 4);

    freeMemArr(table, m);

    return 0;
}
