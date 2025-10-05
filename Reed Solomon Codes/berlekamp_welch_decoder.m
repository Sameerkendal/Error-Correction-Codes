function decoded_msg = berlekamp_welch_decoder(evalpts,y,n, k, tc, prim_poly)
    
        degQ = (k-1)+tc;
        degE = tc;

        NumUnknowns = degQ +degE+1;
        A   = gf(zeros(n, degQ+degE+1), 8, prim_poly);  
        B = gf(zeros(n, 1), 8, prim_poly);

        %% MATRIX FORMATION
        for row  =1:n
            for col =1:degQ+1
                A(row,col) = evalpts(row)^(col-1);
            end
            for col = degQ+2:NumUnknowns
                A(row,col)  = y(row)*(evalpts(row)^(col-(degQ+2)));
            end
        end
        
        for row = 1:n
            B(row)= y(row)*(evalpts(row)^tc);
        end

        %% Gaussian Elimination

        Aug = [A B];

        pivot_cols = [];
        row = 1;  
        
        for col = 1:NumUnknowns
            
            col_slice = Aug(row:n, col);
            col_x = col_slice.x;                
            nonzero_rel = find(col_x ~= 0, 1);   
        
            if isempty(nonzero_rel)
                continue;
            end
        
            pivot_row = nonzero_rel + row - 1;
        
           
            if pivot_row ~= row
                temp = Aug(row, :);
                Aug(row, :) = Aug(pivot_row, :);
                Aug(pivot_row, :) = temp;
            end
        
           
            Aug(row, :) = Aug(row, :) / Aug(row, col);
        
            for r = 1:n
                if r ~= row
                    Aug(r, :) = Aug(r, :) - Aug(r, col) * Aug(row, :);
                end
            end
        
            pivot_cols = [pivot_cols col];
            row = row + 1;
        
            if row > n
                break;
            end
        end
        
      
        free_cols = setdiff(1:NumUnknowns, pivot_cols);
       
        x = gf(zeros(NumUnknowns,1), 8, prim_poly);
        for i = 1:length(pivot_cols)
            x(pivot_cols(i)) = Aug(i, NumUnknowns+1);
        end
        
        if ~isempty(free_cols)
            disp('Infinite solutions: Free variables present.');
        end
        

        %% Coefficients of Q and E 
        x = Aug(:, end);                    

        q = x(1:degQ+1).';                   
        e = x(degQ+2:degQ+1+degE).';                 
        Q = q;                              
        E = [e, gf(1, 8, prim_poly)];

        %% Divison
        if all(Q.x==0)
            Message = gf(zeros(1,100),8,prim_poly);
        else
            [Message, rem] = gf_poly_div(Q, E);
        end

        Message = trim_high(Message);
        if length(Message) < k
            Message = [Message, gf(zeros(1, k - length(Message)), 8, prim_poly)];
        elseif length(Message) > k
            Message = Message(1:k);
        end

        %%
        G = gf(zeros(k, n), 8);

        for ri=1:k
            for ci=1:n
                G(ri, ci)=gf(ci-1, 8)^(ri-1);
            end
        end
        decoded_msg = Message * G(1:k,1:k);
end

function [q, r] = gf_poly_div(num, den)
    num = trim_high(num);
    den = trim_high(den);
    if all(den == 0)
        error('Division by zero polynomial'); 
    end

    m = num.m; 
    pp = num.prim_poly;

    Nh = fliplr(num);
    Dh = fliplr(den);

    degN = length(Nh) - 1;
    degD = length(Dh) - 1;

    if degN < degD
        q = gf(zeros(1,degN), m, pp);
        r = num;
        return;
    end

    Qh = gf(zeros(1, degN - degD + 1), m, pp);
    R  = Nh;

    for t = 1:(degN - degD + 1)
        i = t;                               
        if R(i) ~= 0
            Qh(t) = R(i) / Dh(1);             
            for j = 1:length(Dh)
                R(i + j - 1) = R(i + j - 1) - Qh(t) * Dh(j);
            end
        else
            Qh(t) = gf(0, m, pp);
        end
    end

    q = fliplr(Qh);                           
    r = fliplr(R);
    r = trim_high(r);
end



function p = trim_high(p)
    while length(p) > 1 && p(end) == 0
        p = p(1:end-1);
    end
end
