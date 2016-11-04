function hash = BuildHash(Key,Value,KeyString,ValueString)

% This function is to build hash in Matlab.
% Key or Value should be array or cell array.
% Input:Key,Value,KeyString(0/1),ValueString(0/1);
% Note£º0 means NUMBERL;1means STRING.
% Songpeng Zu
% 120616

%% Main Function
hash = java.util.Hashtable;
klength = length(Key);
vlength = length(Value);
if nargin < 3
    KeyString = 0;
end
if nargin < 4
    ValueString = 0;
end

if klength ~= vlength
    error('Key and Value lengths are not equal.');
else
    if KeyString==0 
        if ValueString == 0
            for i = 1:klength
                hash.put(Key(i),Value(i));
            end
        elseif ValueString == 1
            for i =1:klength
                hash.put(Key(i),Value{i});
            end
        else
            error('ValueString should be 1 or 0.');
        end
    elseif KeyString == 1
        if ValueString == 0
            for i = 1:klength
                hash.put(Key{i},Value(i));
            end
        elseif ValueString ==1
            for i =1:klength
                hash.put(Key{i},Value{i});
            end
        else
            error('ValueString should be 1 or 0.');
        end
    else
        error('KeyString should be 1 or 0.');
    end
end

end

