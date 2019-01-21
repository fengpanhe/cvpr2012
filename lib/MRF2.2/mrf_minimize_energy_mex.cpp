/*
 * Filename: mrf_minimize_energy_mex.cpp
 * Project: MRF2.2
 * Created Date: Tuesday January 15th 2019
 * Author: Feng Panhe
 * -----
 *  Last Modified:
 *  Modified By:
 * -----
 * Copyright (c) 2019 Feng Panhe
 */

#include "mex.hpp"
#include "mexAdapter.hpp"

using matlab::mex::ArgumentList;
using namespace matlab::data;

#include "TRW-S.h"
#include "mrf.h"

struct ORDER_LABEL {
    int abc = 0;
    int acb = 1;
    int bac = 2;
    int bca = 3;
    int cab = 4;
    int cba = 5;
} LABEL;

struct TJunction {
    double sorces[3];
    int eids[3];
};

MRF::CostVal* D;
MRF::CostVal* V;

double** infoMatrix = NULL;
struct TJunction* tjs = NULL;
int tid_index = 0;
int edgeid_index = 1;
int score_index = 2;
double penalties = 0;

MRF::CostVal dCost(int pix, int label)
{
    MRF::CostVal cost = 0;
    int index = pix * 3;
    double a = infoMatrix[index][score_index];
    double b = infoMatrix[index + 1][score_index];
    double c = infoMatrix[index + 2][score_index];
    switch (label) {
    case LABEL.abc:
        cost = a * 2 - c * 2;
        break;
    case LABEL.acb:
        cost = a * 2 - b * 2;
        break;
    case LABEL.bac:
        cost = b * 2 - c * 2;
        break;
    case LABEL.bca:
        cost = b * 2 - a * 2;
        break;
    case LABEL.cab:
        cost = c * 2 - b * 2;
        break;
    case LABEL.cba:
        cost = c * 2 - a * 2;
        break;

    default:
        break;
    }
    return cost / 2;
}

MRF::CostVal fnCost(int pix1, int pix2, int label1, int label2)
{

    MRF::CostVal answer = 0;

    int row_indexs1[3];
    int row_indexs2[3];
    row_indexs1[0] = pix1 * 3;
    row_indexs1[1] = row_indexs1[0] + 1;
    row_indexs1[2] = row_indexs1[0] + 2;

    row_indexs2[0] = pix2 * 3;
    row_indexs2[1] = row_indexs2[0] + 1;
    row_indexs2[2] = row_indexs2[0] + 2;

    int equal_edgeid_index[3];

    for (int i = 0; i < 3, i++) {
        equal_edgeid_index[i] = -1;
        for (int j = 0; j < 3; j++) {
            if (infoMatrix[row_indexs1[i]][edgeid_index] == infoMatrix[row_indexs2[j]][edgeid_index]) {
                equal_edgeid_index[i] = j;
                break;
            }
        }
    }

    for (int i = 0; i < 3; i++) {
        if (equal_edgeid_index[i] != -1) {
            bool flag1 = false;
            bool flag2 = false;

            switch (i) {
            case 0:
                if (label1 == LABEL.abc || label1 == LABEL.acb || label1 == LABEL.cab) {
                    flag1 = true;
                }
                break;
            case 1:
                if (label1 == LABEL.abc || label1 == LABEL.bac || label1 == LABEL.bca) {
                    flag1 = true;
                }
                break;
            case 2:
                if (label1 == LABEL.bca || label1 == LABEL.cba || label1 == LABEL.cab) {
                    flag1 = true;
                }
                break;

            default:
                break;
            }

            switch (equal_edgeid_index[i]) {
            case 0:
                if (label2 == LABEL.abc || label2 == LABEL.acb || label2 == LABEL.cab) {
                    flag2 = true;
                }
                break;
            case 1:
                if (label2 == LABEL.abc || label2 == LABEL.bac || label2 == LABEL.bca) {
                    flag2 = true;
                }
                break;
            case 2:
                if (label2 == LABEL.bca || label2 == LABEL.cba || label2 == LABEL.cab) {
                    flag2 = true;
                }
                break;

            default:
                break;
            }

            if (flag1 == flag2) {
                answer += penalties;
            }
        }
    }

    return answer;
}

EnergyFunction* generate_DataFUNCTION_SmoothGENERAL_FUNCTION()
{
    DataCost* data = new DataCost(dCost);
    SmoothnessCost* smooth = new SmoothnessCost(fnCost);
    EnergyFunction* energy = new EnergyFunction(data, smooth);

    return energy;
}

class MexFunction : public matlab::mex::Function {
    void operator()(matlab::mex::ArgumentList outputs,
        matlab::mex::ArgumentList inputs)
    {
        checkArguments(outputs, inputs);

        penalties = inputs[0][0];
        matlab::data::TypedArray<double> inMatrix = std::move(inputs[1]);

        int row_num = inputs[2][0];
        int col_num = inputs[2][1];

        int t_junction_num = row_num / 3;
        tjs = new TJunction[t_junction_num];

        for (int i = 0; i < t_junction_num; i++) {
            int index0 = i * 3;
            int index1 = index0 + 1;
            int index2 = index0 + 2;
            if (inMatrix[index0][1] != inMatrix[index1][1] || inMatrix[index0][1] != inMatrix[index2][1]) {
                matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
                    0, std::vector<matlab::data::Array>({ factory.createScalar("Tid mismatch") }));
            }
            if (inMatrix[index0][0] != inMatrix[index1][0] || inMatrix[index0][0] != inMatrix[index2][0]) {
                matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
                    0, std::vector<matlab::data::Array>({ factory.createScalar("Imageid mismatch") }));
            }
            tjs->eids[0] = inMatrix[index0][2];
            tjs->eids[1] = inMatrix[index1][2];
            tjs->eids[2] = inMatrix[index2][2];
            tjs->sorces[0] = inMatrix[index0][3];
            tjs->sorces[1] = inMatrix[index1][3];
            tjs->sorces[2] = inMatrix[index2][3];
        }

        int imid = inMatrix[0][0];
        int im_start_index = 0;
        int t_j_count = 0;
        for (int i = 0; i < row_num; i++) {
            if (inMatrix[i][0] != imid) {
                this->minimize_use_TRWS();
            } else {
                t_j_count++;
            }
        }

        // infoMatrix = new double*[row_num + 1];
        // for (int i = 0; i < row_num; i++) {
        //     infoMatrix[i] = new double[col_num + 1];
        // }

        // for (int i = 0; i < row_num; i++) {
        //     for (int j = 0; j < col_num; j++) {
        //         infoMatrix[i][j] = inMatrix[i][j];
        //     }
        // }
        // this->minimize_use_TRWS();
        // for (int i = 0; i < row_num; i++) {
        //     delete[] infoMatrix[i];
        // }
        // delete[] infoMatrix;
    }

    void checkArguments(matlab::mex::ArgumentList outputs,
        matlab::mex::ArgumentList inputs)
    {
        // std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        // matlab::data::ArrayFactory factory;

        if (inputs.size() != 3) {
            matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
                0, std::vector<matlab::data::Array>({ factory.createScalar("Two inputs required") }));
        }

        if (inputs[0].getType() != matlab::data::ArrayType::DOUBLE || inputs[0].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE || inputs[0].getNumberOfElements() != 1) {
            matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
                0, std::vector<matlab::data::Array>({ factory.createScalar("Input penalties must be a scalar") }));
        }

        if (inputs[1].getType() != matlab::data::ArrayType::DOUBLE || inputs[1].getType() == matlab::data::ArrayType::COMPLEX_DOUBLE) {
            matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
                0, std::vector<matlab::data::Array>({ factory.createScalar("Input matrix must be type double") }));
        }

        if (inputs[1].getDimensions().size() != 2) {
            matlabPtr->feval(matlab::engine::convertUTF8StringToUTF16String("error"),
                0, std::vector<matlab::data::Array>({ factory.createScalar("Input must be m-by-n dimension") }));
        }
    }

    void displayOnMATLAB(std::ostringstream& stream)
    {
        // Pass stream content to MATLAB fprintf function
        matlabPtr->feval(u"fprintf", 0,
            std::vector<Array>({ factory.createScalar(stream.str()) }));
        // Clear stream buffer
        stream.str("");
    }

    void minimize_use_TRWS(sizeX, sizeY, numLabels, energy)
    {
        MRF* mrf;
        EnergyFunction* energy;
        MRF::EnergyVal E;
        double lowerBound;
        float t, tot_t;
        int iter;

        energy = generate_DataFUNCTION_SmoothGENERAL_FUNCTION();

        // printf("\n*******Started TRW-S *****\n");
        stream << "\n*******Started TRW-S *****\n";
        this.displayOnMATLAB(stream);

        mrf = new TRWS(sizeX, sizeY, numLabels, energy);

        // can disable caching of values of general smoothness function:
        // mrf->dontCacheSmoothnessCosts();

        mrf->initialize();
        mrf->clearAnswer();

        E = mrf->totalEnergy();
        printf("Energy at the Start= %g (%g,%g)\n", (float)E,
            (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

        tot_t = 0;
        for (iter = 0; iter < 10; iter++) {
            mrf->optimize(10, t);

            E = mrf->totalEnergy();
            lowerBound = mrf->lowerBound();
            tot_t = tot_t + t;
            printf("energy = %g, lower bound = %f (%f secs)\n", (float)E,
                lowerBound, tot_t);
        }
        for (int pix = 0; pix < sizeX * sizeY; pix++)
            printf("Label of pixel %d is %d", pix, mrf->getLabel(pix));

        delete mrf;
    }

private:
    // Pointer to MATLAB engine to call fprintf
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();

    // Factory to create MATLAB data arrays
    matlab::data::ArrayFactory factory;

    // Create an output stream
    std::ostringstream stream;
};