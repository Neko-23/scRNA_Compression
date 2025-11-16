import matplotlib.pyplot as plt
import numpy as np
import os

def create_bar_charts(path, sample, k):
    # Code retrieved from https://stackoverflow.com/a/1392549
    def get_size(start_path = '.'):
        total_size = 0
        for dirpath, dirnames, filenames in os.walk(start_path):
            for f in filenames:
                fp = os.path.join(dirpath, f)
                # skip if it is symbolic link
                if not os.path.islink(fp):
                    total_size += os.path.getsize(fp)
        return round(total_size/1000000,2)

    types = zip(['', '_delta_extension', '_count_extension', '_count_and_delta_extension'], ['', ' further delta.csv compression', ' further counts.csv compression', ' further counts.csv and delta.csv compression'])
    for directory_extension, title_extension in types:
        storage_data = []
        storage_data.append((get_size(path + '/high_level_compress' + directory_extension), 'H'+str(k)+' (ours)', 'green'))
        storage_data.append((round(os.path.getsize(path + '/high_level_compress'+ directory_extension + '.tar.gz')/1000000,2), 'H'+str(k)+'.tar.gz (ours)', 'red'))
        storage_data.append((get_size(path + '/low_level_compress' + directory_extension), 'L'+str(k)+' (ours)', 'green'))
        storage_data.append((round(os.path.getsize(path + '/low_level_compress' + directory_extension + '.tar.gz')/1000000,2), 'L'+str(k)+'.tar.gz (ours)', 'red'))
        storage_data.append((round(get_size('data/sample'+str(sample)+'/csc_nozip')), 'CSC', 'blue'))
        storage_data.append((round(os.path.getsize('data/sample'+str(sample)+'/csc.npz')/1000000,2), 'CSC.tar.gz', 'purple'))
        storage_data.append((round(get_size('data/sample'+str(sample)+'/csr_nozip')), 'CSR', 'blue'))
        storage_data.append((round(os.path.getsize('data/sample'+str(sample)+'/csr.npz')/1000000,2), 'CSR.tar.gz', 'purple'))
        storage_data.append((round(os.path.getsize('data/sample'+str(sample)+'/matrix.mtx')/1000000,2), 'MTX', 'blue'))
        storage_data.append((round(os.path.getsize('data/sample'+str(sample)+'/matrix.tar.gz')/1000000,2), 'MTX.tar.gz', 'purple'))

        storage_data.sort(reverse=True)
        categories = [x[1] for x in storage_data]
        values = [x[0] for x in storage_data]
        colors = [x[2] for x in storage_data]

        # Create the bar chart
        plt.figure(figsize=(15, 8))
        plt.bar(categories, values, color = colors)

        # Add values on top of bars
        for i in range(len(categories)):
            plt.text(i, values[i], str(values[i]) + ' MB', ha='center', va='bottom')

        # Add labels and title
        plt.xlabel('Compression Formats', fontsize=15, labelpad=15)
        plt.ylabel('Storage (MB)', fontsize=15, labelpad=15)
        plt.title('Storage of scRNA Matrix '+str(sample)+' over the Type of Compression Used (ours: k='+str(k) + title_extension + ')', fontsize=20, pad=15)

        # Adjust layout to prevent labels from overlapping
        plt.tight_layout()

        # Save the fig
        plt.savefig(path + '/storage_comparisons' + directory_extension + '.png', dpi=300, bbox_inches='tight')

def create_pie_chart(path, sample, k):
    data = []
    for prefix in ['high', 'low']:
        for extension in ['', '_count_extension', '_delta_extension', '_count_and_delta_extension']:
            cur = os.path.join(path, prefix+'_level_compress'+extension)
            sizes = []
            labels = []
            for name in os.listdir(cur):
                sizes.append(round(os.path.getsize(os.path.join(cur, name))/1000000, 2))
                if sizes[-1] == 0:
                    sizes = sizes[:-1]
                else:
                    labels.append(name.replace('_', '\n'))
            title = prefix.capitalize() + ' Level Compression'
            if extension:
                title += '\n'
                tmp = extension.split('_')
                for ele in tmp:
                    if ele:
                        title += ' '
                        if ele != 'and':
                            title += ele.capitalize()
                        else:
                            title += ele
            data.append((sizes, labels, title))
        
    # Create a figure and a 2x4 grid of subplots (axes)
    fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(14, 7))

    # Flatten the axes array for easier iteration
    axes = axes.flatten()

    total_mb = 0
    def autopct_mb_format(pct):
        # Calculate the actual value from the percentage
        actual_value = (pct / 100) * total_mb
        # Format the value to display as an integer with "MB" suffix
        return f'{pct:.1f}%\n({actual_value:.2f} MB)'

    # Iterate through the data and plot each pie chartx
    colors = plt.cm.Set2.colors
    for i, ax in enumerate(axes):
        if i < len(data):  # Ensure we don't go out of bounds if data is less than 8
            total_mb = sum(data[i][0])
            ax.pie(data[i][0], labels=data[i][1], colors=colors, autopct=autopct_mb_format, pctdistance=0.75, startangle=90, textprops={'fontsize': 8})
            ax.set_title(data[i][2])
            ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    
    fig.suptitle('Storage of scRNA Matrix '+str(sample)+' over the Type of Compression Used for k='+str(k)+' Clusters', fontsize=16)
    plt.tight_layout()  # Adjusts subplot parameters for a tight layout
    plt.savefig(os.path.join(path, 'storage_pie_chart.png'))
    plt.show()
