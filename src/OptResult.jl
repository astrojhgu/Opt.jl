module OptResult

struct MaxIterReached end
struct Finished end

ResultType=Union{Finished, MaxIterReached}

end
