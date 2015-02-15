function Value = GetHash2(hash,Key,KeyString,ValueString)

% This function is to build hash in Matlab.
% Key or Value should be array or cell array.
% Input:Hash,Value,KeyString(0/1),ValueString(0/1);
% Note£º0 means NUMBERL;1means STRING.
%       Compare with GetHash, we consider the situation of Key not in the
%       hash.
% Songpeng Zu /*zusongpeng@gmail.com*/
% 140227



%% Main Function
    klength = length(Key);
    if nargin < 3
        KeyString = 0;
    end
    if nargin < 4
        ValueString = 0;
    end
    if KeyString==0 
        if ValueString == 0
            Value = zeros(klength,1);
            for i = 1:klength
                if hash.containsKey(Key(i))
                    Value(i) = hash.get(Key(i));
                else
                    Value(i) = nan;
                end
            end
        elseif ValueString == 1
            Value{klength} = nan;
            for i =1:klength
                if hash.containsKey(Key(i))
                    Value{i} = hash.get(Key(i));
                else
                    Value{i} = nan;
                end
            end
        else
            error('ValueString should be 1 or 0.');
        end
    elseif KeyString == 1
        if ValueString == 0
            Value = zeros(klength,1);
            for i = 1:klength
                if hash.containsKey(Key(i))
                    Value(i) = hash.get(Key(i));
                else
                    Value(i) = nan;
                end
            end
        elseif ValueString ==1
            Value{klength} = nan;
            for i =1:klength
                if hash.containsKey(Key{i})
                    Value{i} = hash.get(Key{i});
                else
                    Value{i} = nan;
                end
            end
        else
            error('ValueString should be 1 or 0.');
        end
    else
        error('KeyString should be 1 or 0.');
    end


end


