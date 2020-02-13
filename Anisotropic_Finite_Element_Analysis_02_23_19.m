%%Finite Element Analysis Model for Loaded Members
%%The elements used in this code can be of varying length, diameter etc.

disp('Choose a model to run for the finite element analysis.');
disp('1: axial loading with both ends fixed');
disp('2: axial loading with one end fixed');
disp('3: torsional loading with one end fixed');
disp('4: torsional loading with both ends fixed');
n = input('input: ');
YN = 0;
switch n
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    case 1 %axial loading with both ends fixed
        disp('you chose the 1st case')
        elements = 5;
        lValues = [5,5,8,5,5];
        dValues = [6,6,6,3,3];
        aValues = pi.*(dValues./2).^2;
        eValues = 10^6.*[30,30,30,10,10];
        fValues = [0,0,-5000,9000,0,0];
        YN = input('Would you like to use the default values? (0=no/1=yes): ');
        if YN==0
            elements = input('number of elements: ');
            fValues = input('forces at each node in row matrix form: ');
            lValues = input('length of each element in row matrix form: ');
            aValues = input('area of each element in row matrix form: ');
            eValues = input('modulus of elasticity of each element in row matrix form: ');            
        end
            uValues = zeros(1,elements+1);
            aMatrix = zeros(elements+1,elements+1);
            cMatrix = zeros(elements-1,elements-1);
            kValues = aValues.*eValues./lValues;
            disp('kValues')
            disp(kValues)
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    case 2 %axial loading with one end fixed
        disp('you chose the 2nd case')
        elements = 5;
        lValues = [5,5,10,6,6];
        aValues = [8,8,6,4,4];
        eValues = 10^6.*[30,30,10,15,15];
        fValues = [0,0,-500,0,0,2000];
        YN = input('Would you like to use the default values? (0=no/1=yes): ');
        if YN==0
            elements = input('number of elements: ');
            fValues = input('forces at each node in row matrix form: ');
            lValues = input('length of each element in row matrix form: ');
            aValues = input('area of each element in row matrix form: ');
            eValues = input('modulus of elasticity of each element in row matrix form: ');            
        end
            uValues = zeros(1,elements+1);
            aMatrix = zeros(elements+1,elements+1);
            cMatrix = zeros(elements,elements);
            kValues = aValues.*eValues./lValues;
            disp('kValues')
            disp(kValues)
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    case 3 %torsional loading with one end fixed  
        disp('you chose the 3rd case')
            elements = 5;
            lValues = [10,10,24,10,10];
            tValues = 10^3.*[0,-12,0,4,0,-10];
            doValues = [6,6,4,2,2];
            diValues = [3,3,0,0,0];
            gValues = 10^6.*[11.5,11.5,3.7,5.6,5.6];
        YN = input('Would you like to use the default values? (0=no/1=yes): ');
        if YN==0
            elements = input('number of elements: ');
            tValues = input('torques at each node in row matrix form: ');
            lValues = input('length of each element in row matrix form: ');
            doValues = input('outer diameter of each element in row matrix form: ');
            diValues = input('inner diameter of each element in row matrix form: ');
            gValues = input('modulus of rigidity of each element in row matrix form: ');            
        end
            uValues = zeros(1,elements+1);
            aMatrix = zeros(elements+1,elements+1);
            cMatrix = zeros(elements,elements);
            jValues = pi.*(doValues.^4-diValues.^4)./32;
            kValues = gValues.*jValues./lValues;
            disp('kValues')
            disp(kValues)
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    case 4 %torsional loading with both ends fixed  
        disp('you chose the 4th case')
            elements = 5;
            lValues = [20,12,12,10,10];
            tValues = 10^3.*[0,-12,0,4,0,0];
            doValues = [6,4,4,2,2];
            diValues = [3,0,0,0,0];
            gValues = 10^6.*[11.5,3.7,3.7,5.6,5.6];
        YN = input('Would you like to use the default values? (0=no/1=yes): ');
        if YN==0
            elements = input('number of elements: ');
            tValues = input('torques at each node in row matrix form: ');
            lValues = input('length of each element in row matrix form: ');
            doValues = input('outer diameter of each element in row matrix form: ');
            diValues = input('inner diameter of each element in row matrix form: ');
            gValues = input('modulus of rigidity of each element in row matrix form: ');            
        end
            uValues = zeros(1,elements+1);
            aMatrix = zeros(elements+1,elements+1);
            cMatrix = zeros(elements-1,elements-1);
            jValues = pi.*(doValues.^4-diValues.^4)./32;
            kValues = gValues.*jValues./lValues;
            disp('kValues')
            disp(kValues)
    otherwise
        disp('invalid entry')
end
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    %putting together the assembly matrix
    for i=1:elements+1
        if i==1
            aMatrix(i,i)=kValues(i);
            aMatrix(i,i+1)=-kValues(i);
            aMatrix(i+1,i)=-kValues(i);
        elseif i==elements+1
            aMatrix(i,i)=kValues(i-1);
        else
            aMatrix(i,i)=kValues(i-1)+kValues(i);
            aMatrix(i,i+1)=-kValues(i);
            aMatrix(i+1,i)=-kValues(i);
        end
    end
    disp('assembly matrix based on global constraints')
    disp(aMatrix)
    
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
if n==1 %axial loading with both ends fixed
    %condensed matrix
    for i=1:elements-1
        for j=1:elements-1
        cMatrix(i,j)=aMatrix(i+1,j+1);
        end
    end
    disp('condensed matrix')
    disp(cMatrix)
    %case 1 - solving for displacements
    C1uValues = zeros(1,elements-1);
    C1fValues = zeros(1,elements-1);
    for i=1:elements-1
        C1fValues(i)=fValues(i+1);
    end
    disp('f values for case 1')
    disp(C1fValues)
    C1uValues = inv(cMatrix)*C1fValues';
    disp('u values for case 1')
    disp(C1uValues)
elseif n==2 %axial loading with one end fixed
    %condensed matrix
    for i=1:elements
        for j=1:elements
        cMatrix(i,j)=aMatrix(i+1,j+1);
        end
    end
    disp('condensed matrix')
    disp(cMatrix)
    %case 2 - solving for displacements
    C2uValues = zeros(1,elements);
    C2fValues = zeros(1,elements);
    for i=1:elements
        C2fValues(i)=fValues(i+1);
    end
    disp('f values for case 2')
    disp(C2fValues)
    C2uValues = inv(cMatrix)*C2fValues';
    disp('u values for case 2')
    disp(C2uValues)
elseif n==3 %torsional loading with one end fixed
    %condensed matrix
    for i=1:elements
        for j=1:elements
        cMatrix(i,j)=aMatrix(i+1,j+1);
        end
    end
    disp('condensed matrix')
    disp(cMatrix)
    %case 3 - solving for angles of twist
    C3tauValues = zeros(1,elements);
    C3tValues = zeros(1,elements);
    for i=1:elements
        C3tValues(i)=tValues(i+1);
    end
    disp('t values for case 3')
    disp(C3tValues)
    C3tauValues = inv(cMatrix)*C3tValues';
    disp('tau values for case 3')
    disp(C3tauValues)
elseif n==4 %torsional loading with both ends fixed
    %condensed matrix
    for i=1:elements-1
        for j=1:elements-1
        cMatrix(i,j)=aMatrix(i+1,j+1);
        end
    end
    disp('condensed matrix')
    disp(cMatrix)
    %case 4 - solving for angles of twist
    C4tauValues = zeros(1,elements-1);
    C4tValues = zeros(1,elements-1);
    for i=1:elements-1
        C4tValues(i)=tValues(i+1);
    end
    disp('t values for case 4')
    disp(C4tValues)
    C4tauValues = inv(cMatrix)*C4tValues';
    disp('tau values for case 4')
    disp(C4tauValues)
else %invalid entry for the original switch statement
    disp('try running the script again');
end        