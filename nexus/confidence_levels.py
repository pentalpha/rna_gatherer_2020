
def load_confidence(intervals_file):
    with open(intervals_file, 'r') as stream:
        lines = [line.rstrip("\n").split(",") for line in stream.readlines()]
        th_lines = lines[2:]
        confidences = [{} for i in range(len(th_lines[0])-1)]
        for cells in th_lines:
            metric = cells[0]
            for i in range(1,len(cells)):
                confidences[i-1][metric] = float(cells[i])
        '''for i in range(len(confidences)):
            print("Confidence "+ str(i)+": " + str(confidences[i]))'''
        return confidences
    return []