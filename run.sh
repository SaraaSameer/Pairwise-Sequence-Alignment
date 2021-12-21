echo "\n--------------------------------------------------------------------"
echo "\nRunning Serial Program"
echo "\n--------------------------------------------------------------------\n"
g++ -o Serial SequenceAlgo_Serial.cpp 
./Serial

echo "\n\n\n--------------------------------------------------------------------"
echo "\nRunning Using Pthreads"
echo "\n--------------------------------------------------------------------\n"
g++ -o Threads SequenceAlgo_Pthreads.cpp -lpthread
./Threads


echo "\n\n\n--------------------------------------------------------------------"
echo "\nRunning Using OpenMP"
echo "\n--------------------------------------------------------------------\n"
g++ -o OMP SequenceAlgo_OpenMP.cpp -fopenmp
./OMP
