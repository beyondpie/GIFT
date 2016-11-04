function Value = GetHash(hash,Key,KeyString,ValueString)

% This function is to build hash in Matlab.
% Key or Value should be array or cell array.
% Input:Hash,Value,KeyString(0/1),ValueString(0/1);
% Note£º0 means NUMBERL;1means STRING.
% Songpeng Zu
% 120616

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
            for i = 1:klength
                Value(i) = hash.get(Key(i));
            end
        elseif ValueString == 1
            for i =1:klength
                Value{i} = hash.get(Key(i));
            end
        else
            error('ValueString should be 1 or 0.');
        end
    elseif KeyString == 1
        if ValueString == 0
            for i = 1:klength
                Value(i) = hash.get(Key{i});
            end
        elseif ValueString ==1
            for i =1:klength
                Value{i} = hash.get(Key{i});
            end
        else
            error('ValueString should be 1 or 0.');
        end
    else
        error('KeyString should be 1 or 0.');
    end


end


