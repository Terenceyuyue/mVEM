function p = presolve_bounds_from_equalities(p)


    
LU = [p.lb p.ub];

p_F_struc = p.F_struc;
n_p_F_struc_cols = size(p_F_struc,2);

if p.K.f >0

    % Find bounds from sum(xi) = 1, xi>0
    for j = 1:p.K.f
        if p_F_struc(j,1)>0
            [row,col,val] = find(p_F_struc(j,:));
            if all(val(2:end) < 0)
                if all(p.lb(col(2:end)-1)>=0)
                    p.ub(col(2:end)-1) = min( p.ub(col(2:end)-1) , val(1)./abs(val(2:end)'));
                end
            end
        end
    end
   
    A = p.F_struc(1:p.K.f,2:end);
    Ap = max(0,A);
    Am = min(0,A);
    
    for j =  find(sum(p.F_struc(1:p.K.f,2:end) | p.F_struc(1:p.K.f,2:end),2)>1)'
        % Simple x == y
        done = 0;
        b = p_F_struc(j,1);
        if b==0
            [row,col,val] = find(p_F_struc(j,:));
            if length(row) == 2
                if val(1) == -val(2)
                    p.lb(col(1)-1) = max(p.lb(col(1)-1),p.lb(col(2)-1));
                    p.lb(col(2)-1) = max(p.lb(col(1)-1),p.lb(col(2)-1));
                    p.ub(col(1)-1) = min(p.ub(col(1)-1),p.ub(col(2)-1));
                    p.ub(col(2)-1) = min(p.ub(col(1)-1),p.ub(col(2)-1));
                    done = 1;
                elseif val(1) == val(2)
                    p.lb(col(1)-1) = max(p.lb(col(1)-1),-p.ub(col(2)-1));
                    p.lb(col(2)-1) = max(-p.ub(col(1)-1),p.lb(col(2)-1));
                    p.ub(col(1)-1) = min(p.ub(col(1)-1),-p.lb(col(2)-1));
                    p.ub(col(2)-1) = min(-p.lb(col(1)-1),p.ub(col(2)-1));
                    done = 1;
                end
            end
        end
        if ~done
          
            a = p_F_struc(j,2:n_p_F_struc_cols);
            ap = Ap(j,:);
            am = Am(j,:);
            
            p_ub = p.ub;
            p_lb = p.lb;

       
%             find_a_pos = find(a>0);
%             find_a_neg = find(a<0);
%             
%             L_pos = p_lb(find_a_pos);
%             U_pos = p_ub(find_a_pos);
%             
%             L_neg = p_lb(find_a_neg);
%             U_neg = p_ub(find_a_neg);
%             
%             L_pos = repmat(L_pos,1,length(L_pos))-diag(L_pos);
%             L_neg = repmat(L_neg,1,length(L_neg))-diag(L_neg);
%             U_pos = repmat(U_pos,1,length(U_pos))-diag(U_pos);
%             U_neg = repmat(U_neg,1,length(U_neg))-diag(U_neg);
%             
%             newlower_pos = (-b - a(find_a_pos)*U_pos - a(find_a_neg)*L_neg)./a(find_a_pos);
%             %newupper_pos = (-b - am*U - ap*L)/ak;
%                     
            %newlower_neg = (-b - am*U - ap*L)/ak;
            %newupper_neg = (-b - ap*U - am*L)/ak;
            
            find_a = find(a);
            if 1%~any(isinf(p_ub(find_a))) &  ~any(isinf(p_lb(find_a)))
                for k = find_a
                    p_ub_k = p_ub(k);
                    p_lb_k = p_lb(k);
                    if (p_ub_k-p_lb_k) > 1e-8
                        L = p_lb;
                        U = p_ub;
                        L(k) = 0;
                        U(k) = 0;
                        ak = a(k);
                        if ak > 0
                            newlower = (-b - ap*U - am*L)/ak;
                            newupper = (-b - am*U - ap*L)/ak;
                        else
                            newlower = (-b - am*U - ap*L)/ak;
                            newupper = (-b - ap*U - am*L)/ak;
                        end
                        p_ub(k) = min(p_ub_k,newupper);
                        p_lb(k) = max(p_lb_k,newlower);
                    end
                end
                p.ub = p_ub;
                p.lb = p_lb;
            end
        end      
    end
end
close = find(abs(p.lb - p.ub) < 1e-12);
p.lb(close) = (p.lb(close)+p.ub(close))/2;
p.ub(close) = p.lb(close);
if ~isequal(LU,[p.lb p.ub])
    p.changedbounds = 1;
end

% 
% 
% function p = presolve_bounds_from_equalities(p)
% LU = [p.lb p.ub];
% if p.K.f >0
%     
%     % Find bounds from sum(xi) = 1, xi>0
%     for j = 1:p.K.f
%         if p.F_struc(j,1)>0
%             [row,col,val] = find(p.F_struc(j,:));
%             if all(val(2:end) < 0)
%                 if all(p.lb(col(2:end)-1)>=0)
%                     p.ub(col(2:end)-1) = min( p.ub(col(2:end)-1) , val(1)./abs(val(2:end)'));
%                 end
%             end
%         end
%     end
%         
%     for j = 1:p.K.f
%         % Simple x == y
%         done = 0;
%         if p.F_struc(j,1)==0
%             [row,col,val] = find(p.F_struc(j,:));
%             if length(row) == 2
%                 if val(1) == -val(2)
%                     p.lb(col(1)-1) = max(p.lb(col(1)-1),p.lb(col(2)-1));
%                     p.lb(col(2)-1) = max(p.lb(col(1)-1),p.lb(col(2)-1));
%                     p.ub(col(1)-1) = min(p.ub(col(1)-1),p.ub(col(2)-1));
%                     p.ub(col(2)-1) = min(p.ub(col(1)-1),p.ub(col(2)-1));
%                     done = 1;
%                 elseif val(1) == val(2)
%                     p.lb(col(1)-1) = max(p.lb(col(1)-1),-p.ub(col(2)-1));
%                     p.lb(col(2)-1) = max(-p.ub(col(1)-1),p.lb(col(2)-1));
%                     p.ub(col(1)-1) = min(p.ub(col(1)-1),-p.lb(col(2)-1));
%                     p.ub(col(2)-1) = min(-p.lb(col(1)-1),p.ub(col(2)-1));
%                     done = 1;
%                 end               
%             end
%         end
%         if ~done
%             b = p.F_struc(j,1);
%             a = p.F_struc(j,2:end);
%             ap = a.*(a>0);
%             am = a.*(a<0);
%             for k = find(a)
%                 L = p.lb;
%                 U = p.ub;
%                 L(k) = 0;
%                 U(k) = 0;
%                 if a(k) > 0 & (p.ub(k)-p.lb(k)) > 1e-8
%                     newlower = (-b - ap*U - am*L)/a(k);
%                     newupper = (-b - am*U - ap*L)/a(k);
%                     p.ub(k) = min(p.ub(k),newupper);
%                     p.lb(k) = max(p.lb(k),newlower);
%                 end
%             end
%             b = -p.F_struc(j,1);
%             a = -p.F_struc(j,2:end);
%             ap = a.*(a>0);
%             am = a.*(a<0);
%             for k = find(a)
%                 L = p.lb;
%                 U = p.ub;
%                 L(k) = 0;
%                 U(k) = 0;
%                 if a(k) > 0 & (p.ub(k)-p.lb(k)) > 1e-8
%                     newlower = (-b - ap*U - am*L)/a(k);
%                     newupper = (-b - am*U - ap*L)/a(k);
%                     p.ub(k) = min(p.ub(k),newupper);
%                     p.lb(k) = max(p.lb(k),newlower);
%                 end
%             end    
%         end      
%     end
% end
% close = find(abs(p.lb - p.ub) < 1e-12);
% p.lb(close) = (p.lb(close)+p.ub(close))/2;
% p.ub(close) = p.lb(close);
% if ~isequal(LU,[p.lb p.ub])
%     p.changedbounds = 1;
% end