# --------------------------------------------------------------------
#    CAPTURE CHILD PROCESS
with open("i_am_working.txt", "w") as f:
    print("This process has started!", file=f, flush=True)

    import multiprocessing.connection
    try:    import dill as pickle
    except: import pickle

    print("Getting authorization..", file=f, flush=True)
    authorization = bytes('{authorization}','ASCII')

    print("Getting sender at port {port_out}..", file=f, flush=True)
    # Establish the connection receiving data.
    sender = multiprocessing.connection.Client(('localhost', {port_out}), authkey=authorization)

    print("Getting listener at port {port_in}..", file=f, flush=True)
    listener = multiprocessing.connection.Listener(('localhost', {port_in}), authkey=authorization)
    print("accepting connection..", file=f, flush=True)
    listener = listener.accept()

    # Recieve the function, globals, and args, kwargs for the function.
    print("Receiving  globals..", file=f, flush=True)
    globals().update(pickle.loads(listener.recv_bytes()))
    print("Receiving function..", file=f, flush=True)
    function = pickle.loads(listener.recv_bytes())
    print("Receiving 'args' and 'kwargs'..", file=f, flush=True)
    args, kwargs = pickle.loads(listener.recv_bytes())
    print("Running the function..", file=f, flush=True)
    # Run the program safely, send the output.
    try:                     output = function(*args, **kwargs)
    except Exception as exc: output = exc
    print(output, file=f, flush=True)
    # Establish the connection sending output data back.
    sender.send_bytes(pickle.dumps(output))
# --------------------------------------------------------------------
