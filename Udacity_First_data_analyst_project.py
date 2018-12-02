import pandas as pd
from tkinter import filedialog
from matplotlib import pyplot as plt

def plot_data():
    filenames = filedialog.askopenfilenames()
    data_frames = []
    for filename in filenames:
        data = pd.read_csv(filename).ffill()
        data['smooth_temp'] = data['avg_temp'].rolling(5).mean()
        data_frames.append(data)
    city_name = []
    for i in range(1, len(data_frames)):
        city_name.append(data_frames[i].iloc[0]['city'])
    city_name.insert(0, "global")
    if len(city_name) != 6:
        print('This program is for global and 4 cities data file')
    else:
        y_g = data_frames[0]['year']
        temp_g = data_frames[0]['smooth_temp']
        y_c1 = data_frames[1]['year']
        temp_c1 = data_frames[1]['smooth_temp']
        y_c2 = data_frames[2]['year']
        temp_c2 = data_frames[2]['smooth_temp']
        y_c3 = data_frames[3]['year']
        temp_c3 = data_frames[3]['smooth_temp']
        y_c4 = data_frames[4]['year']
        temp_c4 = data_frames[4]['smooth_temp']

        fig, ax = plt.subplots()
        ax.plot(y_g, temp_g, label='Global_Temp')
        ax.plot(y_c1, temp_c1, label=city_name[1])
        ax.plot(y_c2, temp_c2, label=city_name[2])
        ax.plot(y_c3, temp_c3, label=city_name[3])
        ax.plot(y_c4, temp_c4, label=city_name[4])
        ax.legend(loc='upper right')
        ax.set_ylim([0, 35])
        plt.xlabel('Year')
        plt.ylabel('Temperature($^\circ$C)')
        plt.show()
    print(city_name)
    print(data_frames)

plot_data()
