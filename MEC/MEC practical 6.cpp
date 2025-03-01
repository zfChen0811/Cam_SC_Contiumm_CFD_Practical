#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

void read_file(std::vector<std::vector<std::vector<double>>> &data){
    std::ifstream file("/Users/chenzefeng/Desktop/Study/CamSC/MEC/Plasma19 EoS.txt"); // 打开文件
    if (!file) {
        std::cerr << "Failed to open file!" << std::endl;
    }
    
    std::vector<std::vector<double>> currentBlock;      // 当前块的二维表格
    std::string line;
    
    
    while (std::getline(file, line)) {
        // 如果是空行
        if (line.empty()) {
            if (!currentBlock.empty()) {
                data.push_back(currentBlock); // 保存当前块
                currentBlock.clear();         // 清空以准备下一块
            }
            continue;
        }
        
        // 解析非空行
        std::istringstream iss(line);
        std::vector<double> row;
        double value;
        
        while (iss >> value) {
            row.push_back(value); // 将每行数据解析为一个向量
        }
        
        currentBlock.push_back(row); // 将当前行加入当前块
    }
    
    file.close(); // 关闭文件
}

void StoreData(std::vector<std::vector<std::vector<double>>> &data, std::vector<double> &NDNPNS,std::vector<double> &Gamma,std::vector<double> &Heats, std::vector<double> &VibraTemperature, std::vector<double> &Density,std::vector<double> &Pressure,     std::vector<std::vector<double>> &SoundSpeed,std::vector<std::vector<double>> &InterEng, std::vector<std::vector<double>> &Temperature, std::vector<std::vector<double>> &EleCond, std::vector<std::vector<double>> &MassFra, std::vector<std::vector<double>> &ThermCond){
    
    for (size_t i = 0; i < data.size(); ++i) {
        for (const auto& row : data[i]) {
            for (double val : row) {
                if(i == 0){
                    NDNPNS.push_back(val);
                }else if(i == 1){
                    Gamma.push_back(val);
                }else if(i == 2){
                    Heats.push_back(val);
                }else if(i == 3){
                    VibraTemperature.push_back(val);
                }else if(i == 4){
                    Density.push_back(val);
                }else if(i == 5){
                    Pressure.push_back(val);
                }else if(i == 6){
                    SoundSpeed.push_back(row);
                    break;
                }else if(i == 7){
                    InterEng.push_back(row);
                    break;
                }else if(i == 8){
                    Temperature.push_back(row);
                    break;
                }else if(i == 9){
                    EleCond.push_back(row);
                    break;
                }else if(i == 10){
                    MassFra.push_back(row);
                    break;
                }else if(i == 11){
                    ThermCond.push_back(row);
                    break;
                }
            }
        }
    }
}

int main(int argc, char *argv[]){
    int ND, NP, NS;
    int row_nums = 0;
    std::vector<double> NDNPNS;
    std::vector<double> Gamma;
    std::vector<double> Heats;
    std::vector<double> VibraTemperature;
    std::vector<double> Density;
    std::vector<double> Pressure;
    std::vector<std::vector<double>> SoundSpeed;
    std::vector<std::vector<double>> InterEng;
    std::vector<std::vector<double>> Temperature;
    std::vector<std::vector<double>> EleCond;
    std::vector<std::vector<double>> MassFra;
    std::vector<std::vector<double>> ThermCond;
    std::vector<std::vector<std::vector<double>>> data; 
    
    read_file(data);
    StoreData(data, NDNPNS, Gamma, Heats, VibraTemperature, Density, Pressure, SoundSpeed, InterEng, Temperature, EleCond, MassFra, ThermCond);
    
    ND = NDNPNS[0];
    NP = NDNPNS[1];
    NS = NDNPNS[2];
    
    for(int i = 0; i < NP; i++){
        for(int j = 0; j < ND; j++){
            std::cout << SoundSpeed[i][j] << " ";
        }
        std::cout << endl;
    }
}
        