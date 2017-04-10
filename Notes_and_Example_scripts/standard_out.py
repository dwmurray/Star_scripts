def pass_to_std_out(output_text) :
	temp = sys.stdout #store original stdout object for later
	if args.first:
		txt_filename = '{0}/output_prints_{1}_{2}_{3}.txt'.format(output_location, args.start, args.end, 'first')
	elif args.second:
		txt_filename = '{0}/output_prints_{1}_{2}_{3}.txt'.format(output_location, args.start, args.end, 'second')
	elif withBackTracking:
		txt_filename = '{0}/output_prints_{1}_{2}_{3}.txt'.format(output_location, args.start, args.end, 'backwards')
	else:
		txt_filename = '{0}/output_prints_{1}_{2}.txt'.format(output_location, args.start, args.end)
        sys.stdout = open(txt_filename,'a') #redirect all prints to this log file
	print(output_text) #nothing appears at interactive prompt
	sys.stdout.close() #ordinary file object
	sys.stdout = temp #restore print commands to interactive prompt

