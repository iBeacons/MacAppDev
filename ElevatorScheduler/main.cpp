//
//  main.cpp
//  ElevatorScheduler
//
//  Created by Xia Mingjie on 13-11-2.
//  Copyright (c) 2013年 iBeacons.github.com. All rights reserved.
//
//  1、遍历楼层, 时间复杂度为O(N*N)
//
//      目标函数1: min T(i), 所有人走的楼层总数最小值
//      目标函数2: min E(i), 所有人消耗能量最小值
//
//      p[i]: 去i层的人数。
//      |j-i|: 电梯停在第i层时去第j层的人要走的楼层。
//
//      所有人走的楼层总数: T(i)=∑(p[j]*|j-i|), j=1, ..., N
//
//      坐电梯上楼消耗能量: E1(i)=k3*(i-1)*∑p[j], j=1, ..., N
//      走路上下楼消耗能量: E2(i)=∑(p[j]*|j-i|*((j-i)?k1:k2)), j=1, ..., N
//      总的消耗能量: E(i)=E1(i)+E2(i)
//
//
//  2、动态规划，时间复杂度为O(N-1)
//
//      i, j: 1~N
//      p[i]: 去i层的人数。
//      |j-i|: 电梯停在第i层时去第j层的人要走的楼层。
//
//      归纳假设, 假设电梯停在第i层时，
//                  所有人走的楼层总数为Q(i), 消耗能量总数E(i);
//                  要到1至(i-1)层的人数总数为S1, 要到i层的人数总数为S2, 要到(i+1)至N层的人数总数为S3;
//                  /// 对应的人群消耗能量分别为W1, W2, W3.(注：Q(1)= ∑(p[j]*|j-1|), j=1, ..., N;
//                  /// E(i)=W1+W2+W3)
//
//              那么电梯停在第(i+1)层时, 要到1至i层的人需要多下一层楼，他们走的楼层总数增加S1+S2, 消耗能量增加
//                  (S1+S2)*(k2+k3)， 要到(i+1)至N层的人走可
//                  以少上一层楼， 他们走的楼层总数减少S3, 消耗能量变化
//                  S3*(k3-k1),(实际上考虑k3<k1的情况，消耗能量减少S3*(k1-k3)）
//
//                  此时所有人走的楼层总数Q(i+1)=Q(i)+S1+S2-S3;
//                  消耗能量总数E(i+1)=E(i)+(S1+S2)*(k2+k3)－S3*(k1-k3)；
//
//       可以看出随着电梯所停楼层i的增加, S1+S2也是不断增大的，而S3不断减小的。
//                      对应的(S1+S2)*(k2+k3)不断增大，考虑k3<k1的情况，S3*(k1-k3)不断减小。
//
//       S1=∑p[k], k=1，..., i-1
//       S2=p[i]
//       S3=∑p[k], k=(i+1), ..., N
//       // 电梯停在第1层时，S1=0, S2=p[1], S3=∑p[k], k=2~N
//       当S1+S2<S3时， Q(i)是递减的；
//       当S1+S2>S3时， Q(i)是递增的；
//
//       所有人走的楼层总数的最小值在S1+S2<S3的最后一次取得；
//
//
//       ////   W1=k3*∑(p[k]*(i-1))+k2*∑p[k]*(i-k), k=1, ..., i-1
//       ////   W2=k3*p[i]*(i-1)
//       ////   w3=k3*∑(p[k]*(i-1))+k1*∑p[k]*(k-i), k=i+1, ..., N
//       // 电梯停在第1层时，W1=0, W2=0, W3=k1*∑p[k]*(k-1), k=2~N
//
//       当(S1+S2)*(k2+k3)<S3*(k1-k3)时， E(i)是递减的；
//       当(S1+S2)*(k2+k3)>S3*(k1-k3)时， E(i)是递增的；
//
//          所有人消耗能量的最小值在(S1+S2)*(k2+k3)<S3*(k1-k3)的最后一次取得。
//
//
//  3、最佳一致逼近，零点附近给出极值. 时间复杂度为O((N-1)*(1~N);
//
//

#include <iostream>
#include <float.h>

void caculateFloor(const int, const int []);
void caculateEnergy(const int, const int [], const float, const float, const float);

void caculateFloorWithDivision(const int, const int []);
void caculateEnergyWithDivision(const int, const int [], const float, const float, const float);

void weightedBestApproximate(const int, const int []);

int main(int argc, const char * argv[])
{
    int N;
    
    std::cout << "Input number of floors: " << std::endl;
    std::cin >> N;
    
    int p[N+1];
    int n = 1;
    
    std::cout << "Input amount of persons to each floor: " << std::endl;
    while (n <= N) {
        std::cin >> p[n];
        ++n;
    }
    p[0] = 0;
    
    char ch;
    
    std::cout << "Input choice number, 1 for min floors and 2 for min energy." << std::endl;
    std::cin >> ch;
    
    //输入1, 开始计算寻找目标楼层1和楼层总数
    if (ch == '1') {
        caculateFloor(N, p);
        caculateFloorWithDivision(N, p);
        weightedBestApproximate(N, p);
    }
    else if (ch == '2') {   //输入2, 开始计算寻找目标楼层1和消耗能量
        float  k1, k2, k3;
        
        std::cout << "Input energy value k1, k2 and k3: " << std::endl;
        std::cin >> k1 >> k2 >> k3;
        
        caculateEnergy(N, p, k1, k2, k3);
        caculateEnergyWithDivision(N, p, k1, k2, k3);
    }
    else
        std::cout << "input choice can't be accessed" << std::endl;
    
    return 0;
}

void caculateFloor(const int N, const int p[])
{
    int targetFloor = 1;
    int totalClimbStairs = 0;
    long long int minTotalClimbStairs = LONG_LONG_MAX;
    
    int count = 0;
    
    for (int i = 1; i <= N; i++) {
        totalClimbStairs = 0;
        
        for (int j = 1; j <= N; j++) {
            totalClimbStairs += p[j] * abs(j-i);
            ++count;
        }
        
        if (minTotalClimbStairs > totalClimbStairs) {
            minTotalClimbStairs = totalClimbStairs;
            targetFloor = i;
        }
    }
    
    std::cout << "Target floor for climbing stairs: " << targetFloor << std::endl;
    std::cout << "Total number of climbing stairs: " << minTotalClimbStairs << std::endl;
    std::cout << "count: " << count << std::endl;
}

void caculateEnergy(const int N, const int p[], const float k1, const float k2, const float k3)
{
    int targetFloor = 1;
    
    long double totalEnergy = 0.0f;
    long double totalEnergyWithElevator = 0.0f;
    long double totalEnergyWithFeet = 0.0f;
    long double minTotalEnergy = LONG_LONG_MAX*1.0f;
    
    for (int i = 1; i <= N; i++) {
        totalEnergy = 0.0f;
        totalEnergyWithElevator = 0.0f;
        totalEnergyWithFeet = 0.0f;
        
        for (int j = 1; j <= N; j++) {
            totalEnergyWithElevator += k3 * (i-1) * p[j];
            totalEnergyWithFeet +=  p[j] * abs(j-i) * (j>i?k1:k2);
        }
        totalEnergy = totalEnergyWithElevator + totalEnergyWithFeet;
        
        if ((totalEnergy - minTotalEnergy) < DBL_EPSILON) {
            minTotalEnergy = totalEnergy;
            targetFloor = i;
        }
    }
    
    std::cout << "Target floor for energy: " << targetFloor << std::endl;
    std::cout << "Total energy: " << minTotalEnergy << std::endl;
}

void caculateFloorWithDivision(const int N, const int p[])
{
    int s1 = 0;
    int s2 = p[1];
    int s3 = 0;
    
    int targetFloor = 1;
    long long int totalClimbStairs = 0ll;
    
    int count = 0;
    
    for (int k = 2; k <= N; k++) {
        s3 += p[k];
        totalClimbStairs += p[k] * (k-1);
    }
    
    for (int i = 2; i <= N; i++) {
        if (s1+s2 < s3) {
            targetFloor = i;
            totalClimbStairs += s1 + s2 - s3;
            
            s1 += s2;
            s2 = p[i];
            s3 -= p[i];
            
            ++count;
        }
        else
            break;
    }
    
    std::cout << "methodWithDivision-Target floor for climbing stairs: " << targetFloor << std::endl;
    std::cout << "methodWithDivision-Total number of climbing stairs: " << totalClimbStairs << std::endl;
    std::cout << "count: " << count << std::endl;
}

void caculateEnergyWithDivision(const int N, const int p[], const float k1, const float k2, const float k3)
{
    int s1 = 0;
    int s2 = p[1];
    int s3 = 0;
    
    int targetFloor = 1;
    long double totalEnergy = 0.0f;
    
    for (int k = 2; k <= N; k++) {
        s3 += p[k];
        totalEnergy += k1 * (k-1) * p[k];
    }
    
    for (int i = 2; i <= N; i++) {
        if (((s1+s2) * (k2+k3) - s3 * (k1-k3)) < DBL_EPSILON) {
            targetFloor = i;
            totalEnergy += (s1+s2) * (k2+k3) - s3 * (k1-k3);
            
            s1 += s2;
            s2 = p[i];
            s3 -= p[i];
        }
        else
            break;
    }
    
    std::cout << "methodWithDivision-Target floor for energy: " << targetFloor << std::endl;
    std::cout << "methodWithDivision-Total energy: " << totalEnergy << std::endl;
}

void weightedBestApproximate(const int N, const int p[])
{
    int targetFloor = 1;
    int derivative = 0;
    
    int count = 0;
    
    for (int i = 1; i <= N; i++) {
        derivative = 0;
        
        for (int j = 1; j <= N; j++)
            if (j != i) {
                derivative += (i-j)/abs(i-j)*p[j];
                ++count;
            }
        
        if (derivative < 0)
            targetFloor = i;
        else
            break;
    }
    
    long long int totalClimbStairs = 0ll;
    long long int totalClimbStairsLeft = 0ll;
    long long int totalClimbStairsRight = 0ll;
    
    for (int j = 1; j <= N; j++)
    {
        totalClimbStairsLeft += abs(targetFloor-j)*p[j];
        totalClimbStairsRight += abs(targetFloor+1-j)*p[j];
    }
    
    if (totalClimbStairsLeft <= totalClimbStairsRight)
        totalClimbStairs = totalClimbStairsLeft;
    else
    {
        targetFloor = targetFloor + 1;
        totalClimbStairs = totalClimbStairsRight;
    }
    
    std::cout << "methodWithApproximate-Target floor for climbing stairs: " << targetFloor << std::endl;
    std::cout << "methodWithApproximate-Total number for climbing stairs: " << totalClimbStairs << std::endl;
    std::cout << "count: " << count << std::endl;
}

